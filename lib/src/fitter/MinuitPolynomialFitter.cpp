#include <iostream>

#include "D2K3PiError.h"
#include "FitterUtils.h"
#include "MinuitPolynomialFitter.h"
#include "physics.h"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "TF1.h"

PolynomialChiSqFcn::PolynomialChiSqFcn(const std::vector<double>& data,
                                       const std::vector<double>& times,
                                       const std::vector<double>& errors,
                                       const IntegralOptions_t&   integralOptions)
    : MyBaseFcn(data, times, errors, integralOptions)
{
    ;
}

PolynomialChiSqFcn::~PolynomialChiSqFcn()
{
    ;
}

double PolynomialChiSqFcn::operator()(const std::vector<double>& parameters) const
{
    size_t numParams = parameters.size();
    if (numParams != 3) {
        std::cerr << "Require three parameters (a, b, c); got " << numParams << std::endl;
        throw D2K3PiException();
    }

    double chi2 = 0.0;

    if (_integralOptions.integrate) {
        for (size_t i = 0; i < theMeasurements.size(); ++i) {
            double model = Phys::dcsIntegralWithEfficiency(_integralOptions.binLimits[i],
                                                           _integralOptions.binLimits[i + 1],
                                                           parameters,
                                                           _integralOptions.width,
                                                           _integralOptions.efficiencyTimescale,
                                                           _integralOptions.tolerance,
                                                           _integralOptions.maxRefinements) /
                           Phys::cfIntegralWithEfficiency(_integralOptions.binLimits[i],
                                                          _integralOptions.binLimits[i + 1],
                                                          _integralOptions.width,
                                                          _integralOptions.efficiencyTimescale,
                                                          _integralOptions.tolerance,
                                                          _integralOptions.maxRefinements);
            chi2 += std::pow((model - theMeasurements[i]) / theMVariances[i], 2);
        }
    } else {
        for (size_t i = 0; i < theMeasurements.size(); ++i) {
            chi2 += std::pow((Phys::rateRatio(thePositions[i], parameters) - theMeasurements[i]) / theMVariances[i], 2);
        }
    }
    return chi2;
}

MinuitPolynomialFitter::MinuitPolynomialFitter(const FitData_t& fitData, const IntegralOptions_t& integralOptions)
    : MinuitScannerBase(fitData)
{
    _fitFcn =
        std::make_unique<PolynomialChiSqFcn>(_fitData.data, _fitData.binCentres, _fitData.errors, integralOptions);
}

void MinuitPolynomialFitter::fit()
{
    if (!_parameters) {
        std::cerr << "parameters not set; call this->setPolynomialParams before fitting." << std::endl;
        throw D2K3PiException();
    }

    MinuitFitterBase::fit();

    bestFitFunction = std::make_unique<TF1>("Best fit function", "[0] +[1]*x+[2]*x*x");
    bestFitFunction->SetParameters(fitParams.fitParams.data());
}

void MinuitPolynomialFitter::setPolynomialParams(const std::vector<double>& initialParams,
                                                 const std::vector<double>& initialErrors)
{
    _setParams(std::vector<std::string>{"a", "b", "c"}, initialParams, initialErrors);
}
