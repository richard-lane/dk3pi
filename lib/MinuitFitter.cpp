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
        std::cerr << "Require three parameters; got " << numParams << std::endl;
        throw D2K3PiException();
    }

    SimplePolynomialFunction MyPolynomial(parameters[0], parameters[1], parameters[2]);
    double                   chi2 = 0.0;

    for (size_t i = 0; i < theMeasurements.size(); ++i) {
        chi2 += std::pow((MyPolynomial(thePositions[i]) - theMeasurements[i]) / theMVariances[i], 2);
    }
    return chi2;
}
