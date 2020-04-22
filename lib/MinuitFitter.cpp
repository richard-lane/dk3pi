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
#include "MinuitFitter.h"
#include "fitter/MinuitPolynomialFitter.h"
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
        double model =
            Phys::analyticalDcsIntegral(
                _integralOptions.binLimits[i], _integralOptions.binLimits[i + 1], parameters, _integralOptions.width) /
            Phys::analyticalCfIntegral(
                _integralOptions.binLimits[i], _integralOptions.binLimits[i + 1], _integralOptions.width);
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
                                                       const std::vector<double>& errors)
    : BasePolynomialFcn(data, times, errors)
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
        chi2 += std::pow((Phys::rateRatio(thePositions[i], params) - theMeasurements[i]) / theMVariances[i], 2);
    }
    return chi2;
}
