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
#include "util.h"

double fitPolynomial(const DecayParams_t& params, const double time)
{
    return (params.r * params.r) +
           (params.r * (params.y * params.z_re + params.x * params.z_im) * params.width * time) +
           (0.25 * (params.x * params.x + params.y * params.y) * params.width * params.width * time * time);
}

SimplePolynomialFunction::SimplePolynomialFunction(double a, double b, double c) : _a(a), _b(b), _c(c)
{
    ;
}

SimplePolynomialFunction::~SimplePolynomialFunction()
{
    ;
}

double SimplePolynomialFunction::operator()(double x) const
{
    return _a + _b * x + _c * x * x;
}

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
                                       const std::vector<double>& errors)
    : BasePolynomialFcn(data, times, errors)
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

    SimplePolynomialFunction MyPolynomial(parameters[0], parameters[1], parameters[2]);
    double                   chi2 = 0.0;

    for (size_t i = 0; i < theMeasurements.size(); ++i) {
        chi2 += std::pow((MyPolynomial(thePositions[i]) - theMeasurements[i]) / theMVariances[i], 2);
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
        chi2 += std::pow((fitPolynomial(params, thePositions[i]) - theMeasurements[i]) / theMVariances[i], 2);
    }
    return chi2;
}
