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

PolynomialFitFcn::PolynomialFitFcn(const std::vector<double>& meas,
                                   const std::vector<double>& pos,
                                   const std::vector<double>& mvar)
    : _theMeasurements(meas), _thePositions(pos), _theMVariances(mvar), _theErrorDef(1.)
{
    ;
}

PolynomialFitFcn::~PolynomialFitFcn()
{
    ;
}

double PolynomialFitFcn::Up() const
{
    return _theErrorDef;
}

double PolynomialFitFcn::operator()(const std::vector<double>& parameters) const
{
    size_t numParams = parameters.size();
    if (numParams != 3) {
        std::cerr << "Require three parameters; got " << numParams << std::endl;
        throw D2K3PiException();
    }

    SimplePolynomialFunction MyPolynomial(parameters[0], parameters[1], parameters[2]);
    double                   chi2 = 0.0;

    for (size_t i = 0; i < _theMeasurements.size(); ++i) {
        chi2 += (std::pow(MyPolynomial(_thePositions[i]) - _theMeasurements[i], 2)) / _theMVariances[i];
    }
    return chi2;
}

std::vector<double> PolynomialFitFcn::measurements() const
{
    return _theMeasurements;
}

std::vector<double> PolynomialFitFcn::positions() const
{
    return _thePositions;
}

std::vector<double> PolynomialFitFcn::variances() const
{
    return _theMVariances;
}

void PolynomialFitFcn::setErrorDef(double def)
{
    _theErrorDef = def;
}
