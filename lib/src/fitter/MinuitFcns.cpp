/*
 * If we want to perform a fit with Minuit2, we need to define a special class that encapsulates the function that we're
 * using.
 *
 * These are defined here and used by the Fitter.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "D2K3PiError.h"
#include "MinuitFcns.h"
#include "MinuitPolynomialFitter.h"
#include "physics.h"
#include "util.h"

BasePolynomialFcn::BasePolynomialFcn(const std::vector<double>& data,
                                     const std::vector<double>& times,
                                     const std::vector<double>& errors)
    : theMeasurements(data), thePositions(times), theMVariances(errors), theErrorDef(1.)
{
    ;
}

BasePolynomialFcn::~BasePolynomialFcn()
{
    ;
}

double BasePolynomialFcn::Up() const
{
    return theErrorDef;
}

std::vector<double> BasePolynomialFcn::measurements() const
{
    return theMeasurements;
}

std::vector<double> BasePolynomialFcn::positions() const
{
    return thePositions;
}

std::vector<double> BasePolynomialFcn::variances() const
{
    return theMVariances;
}

void BasePolynomialFcn::setErrorDef(double def)
{
    theErrorDef = def;
}

PolynomialChiSqFcn::PolynomialChiSqFcn(const std::vector<double>& data,
                                       const std::vector<double>& times,
                                       const std::vector<double>& errors,
                                       const IntegralOptions_t&   integralOptions)
    : BasePolynomialFcn(data, times, errors), _integralOptions(integralOptions)
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

    for (size_t i = 0; i < theMeasurements.size(); ++i) {
        double model = Phys::dcsRateWithEfficiency(
                           thePositions[i], parameters, _integralOptions.width, 1 / _integralOptions.width) /
                       Phys::cfRateWithEfficiency(thePositions[i], _integralOptions.width, 1 / _integralOptions.width);
        chi2 += std::pow((model - theMeasurements[i]) / theMVariances[i], 2);
    }
    return chi2;
}

PolynomialChiSqFcnNoIntegral::PolynomialChiSqFcnNoIntegral(const std::vector<double>& data,
                                                           const std::vector<double>& times,
                                                           const std::vector<double>& errors)
    : BasePolynomialFcn(data, times, errors)
{
    ;
}

PolynomialChiSqFcnNoIntegral::~PolynomialChiSqFcnNoIntegral()
{
    ;
}

double PolynomialChiSqFcnNoIntegral::operator()(const std::vector<double>& parameters) const
{
    size_t numParams = parameters.size();
    if (numParams != 3) {
        std::cerr << "Require three parameters (a, b, c); got " << numParams << std::endl;
        throw D2K3PiException();
    }

    double chi2 = 0.0;

    for (size_t i = 0; i < theMeasurements.size(); ++i) {
        chi2 += std::pow((Phys::rateRatio(thePositions[i], parameters) - theMeasurements[i]) / theMVariances[i], 2);
    }
    return chi2;
}

DetailedPolynomialChiSqFcn::DetailedPolynomialChiSqFcn(const std::vector<double>& data,
                                                       const std::vector<double>& times,
                                                       const std::vector<double>& errors,
                                                       const IntegralOptions_t&   integralOptions)
    : BasePolynomialFcn(data, times, errors), _integralOptions(integralOptions)
{
    ;
}

DetailedPolynomialChiSqFcn::~DetailedPolynomialChiSqFcn()
{
    ;
}

double DetailedPolynomialChiSqFcn::operator()(const std::vector<double>& parameters) const
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
    for (size_t i = 0; i < theMeasurements.size(); ++i) {
        double model = Phys::dcsRateWithEfficiency(thePositions[i], params, 1 / params.width) /
                       Phys::cfRateWithEfficiency(thePositions[i], params, 1 / params.width);
        chi2 += std::pow((model - theMeasurements[i]) / theMVariances[i], 2);
    }
    return chi2;
}

DetailedChiSqConstrainXYFcn::DetailedChiSqConstrainXYFcn(const std::vector<double>& data,
                                                         const std::vector<double>& times,
                                                         const std::vector<double>& errors,
                                                         const IntegralOptions_t&   integralOptions)
    : BasePolynomialFcn(data, times, errors), _integralOptions(integralOptions)
{
    ;
}

DetailedChiSqConstrainXYFcn::~DetailedChiSqConstrainXYFcn()
{
    ;
}

double DetailedChiSqConstrainXYFcn::operator()(const std::vector<double>& parameters) const
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
    for (size_t i = 0; i < theMeasurements.size(); ++i) {
        double model = Phys::dcsRateWithEfficiency(thePositions[i], params, 1 / params.width) /
                       Phys::cfRateWithEfficiency(thePositions[i], params, 1 / params.width);
        chi2 += std::pow((model - theMeasurements[i]) / theMVariances[i], 2);
    }
    // Introduce constraint by modifying chi squared
    double dx         = parameters[0] - WORLD_AVERAGE_X;
    double dy         = parameters[1] - WORLD_AVERAGE_Y;
    double constraint = (1 / (1 - X_Y_CORRELATION * X_Y_CORRELATION)) *
                        (std::pow(dx / WORLD_AVERAGE_X_ERR, 2) + std::pow(dy / WORLD_AVERAGE_Y_ERR, 2) -
                         2 * X_Y_CORRELATION * dx * dy / (WORLD_AVERAGE_Y_ERR * WORLD_AVERAGE_X_ERR));

    return chi2 + constraint;
}
