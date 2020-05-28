#include <iostream>

#include "D2K3PiError.h"
#include "PhysicalFitter.h"
#include "physics.h"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "TF1.h"

PhysicalFitFcn::PhysicalFitFcn(const std::vector<double>& data,
                               const std::vector<double>& times,
                               const std::vector<double>& errors,
                               const IntegralOptions_t&   integralOptions)
    : MyBaseFcn(data, times, errors, integralOptions)
{
    ;
}

PhysicalFitFcn::~PhysicalFitFcn()
{
    ;
}

double PhysicalFitFcn::operator()(const std::vector<double>& parameters) const
{
    size_t numParams = parameters.size();
    if (numParams != 6) {
        std::cerr << "Require six parameters(x, y, r, z_im, z_re, width); got " << numParams << std::endl;
        throw D2K3PiException();
    }

    DecayParams_t params = DecayParameters{.x     = parameters[0],
                                           .y     = parameters[1],
                                           .r     = parameters[2],
                                           .z_im  = parameters[3],
                                           .z_re  = parameters[4],
                                           .width = parameters[5]};

    double chi2 = 0.0;
    if (_integralOptions.integrate) {
        for (size_t i = 0; i < theMeasurements.size(); ++i) {
            double model = Phys::dcsIntegralWithEfficiency(_integralOptions.binLimits[i],
                                                           _integralOptions.binLimits[i + 1],
                                                           util::expectedParams(params),
                                                           params.width,
                                                           _integralOptions.efficiencyTimescale) /
                           Phys::cfIntegralWithEfficiency(_integralOptions.binLimits[i],
                                                          _integralOptions.binLimits[i + 1],
                                                          params.width,
                                                          _integralOptions.efficiencyTimescale);
            chi2 += std::pow((model - theMeasurements[i]) / theMVariances[i], 2);
        }
    } else {
        for (size_t i = 0; i < theMeasurements.size(); ++i) {
            double model = Phys::rateRatio(thePositions[i], params);
            chi2 += std::pow((model - theMeasurements[i]) / theMVariances[i], 2);
        }
    }
    return chi2;
}

ConstrainXYFcn::ConstrainXYFcn(const std::vector<double>& data,
                               const std::vector<double>& times,
                               const std::vector<double>& errors,
                               const IntegralOptions_t&   integralOptions)
    : MyBaseFcn(data, times, errors, integralOptions)
{
    ;
}

ConstrainXYFcn::~ConstrainXYFcn()
{
    ;
}

double ConstrainXYFcn::operator()(const std::vector<double>& parameters) const
{
    size_t numParams = parameters.size();
    if (numParams != 6) {
        std::cerr << "Require six parameters(x, y, r, z_im, z_re, width); got " << numParams << std::endl;
        throw D2K3PiException();
    }

    DecayParams_t params = DecayParameters{.x     = parameters[0],
                                           .y     = parameters[1],
                                           .r     = parameters[2],
                                           .z_im  = parameters[3],
                                           .z_re  = parameters[4],
                                           .width = parameters[5]};

    std::vector<double> expectedParams = util::expectedParams(params);

    double chi2 = 0.0;
    if (_integralOptions.integrate) {
        for (size_t i = 0; i < theMeasurements.size(); ++i) {
            double model = Phys::dcsIntegralWithEfficiency(_integralOptions.binLimits[i],
                                                           _integralOptions.binLimits[i + 1],
                                                           util::expectedParams(params),
                                                           params.width,
                                                           _integralOptions.efficiencyTimescale) /
                           Phys::cfIntegralWithEfficiency(_integralOptions.binLimits[i],
                                                          _integralOptions.binLimits[i + 1],
                                                          params.width,
                                                          _integralOptions.efficiencyTimescale);
            chi2 += std::pow((model - theMeasurements[i]) / theMVariances[i], 2);
        }
    } else {
        for (size_t i = 0; i < theMeasurements.size(); ++i) {
            chi2 += std::pow((Phys::rateRatio(thePositions[i], params) - theMeasurements[i]) / theMVariances[i], 2);
        }
    }
    // Introduce constraint by modifying chi squared
    double dx         = parameters[0] - WORLD_AVERAGE_X;
    double dy         = parameters[1] - WORLD_AVERAGE_Y;
    double constraint = (1 / (1 - X_Y_CORRELATION * X_Y_CORRELATION)) *
                        (std::pow(dx / WORLD_AVERAGE_X_ERR, 2) + std::pow(dy / WORLD_AVERAGE_Y_ERR, 2) -
                         2 * X_Y_CORRELATION * dx * dy / (WORLD_AVERAGE_Y_ERR * WORLD_AVERAGE_X_ERR));

    return chi2 + constraint;
}

PhysicalFitter::PhysicalFitter(const FitData_t&         fitData,
                               const IntegralOptions_t& integralOptions,
                               const bool               constrainXY)
    : MinuitScannerBase(fitData), _constrainXY(constrainXY)
{
    if (_constrainXY) {
        _fitFcn =
            std::make_unique<ConstrainXYFcn>(_fitData.data, _fitData.binCentres, _fitData.errors, integralOptions);
    } else {
        _fitFcn =
            std::make_unique<PhysicalFitFcn>(_fitData.data, _fitData.binCentres, _fitData.errors, integralOptions);
    }
}

void PhysicalFitter::fit()
{
    // Check that parameters have been set
    if (!_parameters) {
        std::cerr << "Parameters not set" << std::endl;
        throw D2K3PiException();
    }

    // Check we have constrained the right parameters
    if (!_constrainXY) {
        // Check that we have fixed at least three of x, y, Re(Z), Im(Z) and decay width.
        // Copy the parameters we have into a vector + remove r
        std::vector<ROOT::Minuit2::MinuitParameter> paramsExceptR = _parameters->Trafo().Parameters();
        paramsExceptR.erase(paramsExceptR.begin() + 2);

        // Count how many of the relevant parameters are fixed
        short unsigned int numFixedParams{0};
        for (auto it = paramsExceptR.begin(); it != paramsExceptR.end(); ++it) {
            if (it->IsFixed()) {
                ++numFixedParams;
            }
        }

        // Raise if not enough are fixed
        if (numFixedParams < 3) {
            std::cerr << "Must fix at least three of {x, y, width, Im(Z), Re(Z)} for fit to be well defined ("
                      << numFixedParams << " fixed)." << std::endl;
            throw D2K3PiException();
        }
    } else {
        // Check we have fixed at least one of {Re(Z), Im(Z) or width}
        std::vector<ROOT::Minuit2::MinuitParameter> reZ_imZ_width{_parameters->Trafo().Parameters()[3],
                                                                  _parameters->Trafo().Parameters()[4],
                                                                  _parameters->Trafo().Parameters()[5]};
        if (std::all_of(reZ_imZ_width.begin(), reZ_imZ_width.end(), [](ROOT::Minuit2::MinuitParameter& param) {
                return !param.IsFixed();
            })) {
            std::cerr << "Must fix one component of Z or width for fit to be well defined" << std::endl;
        }
    }

    // Use base class implementation to actually perform the fit
    MinuitFitterBase::fit();

    // Create a best-fit function
    bestFitFunction = std::make_unique<TF1>(
        "Best fit function", "[2]*[2] + [2]*([1]*[4] + [0]*[3])*[5]*x + 0.25 * ([0]*[0] + [1]*[1])*[5]*[5]*x*x");
    bestFitFunction->SetParameters(fitParams.fitParams.data());
}

void PhysicalFitter::setFitParams(const std::vector<double>& initialParams, const std::vector<double>& initialErrors)
{
    if (initialParams.size() != 6 || initialErrors.size() != 6) {
        std::cerr << "Must provide 6 initial parameters and errors for x, y, rD, Im(Z), Re(Z) and decay width."
                  << std::endl;
        throw D2K3PiException();
    }
    _setParams(_paramNames, initialParams, initialErrors);
}
