/*
 * Likelihood Fit binned data
 * Intended to be an upgrade of Fitter.h
 */
#ifndef MINUIT_FITTER_H
#define MINUIT_FITTER_H

#include <memory>
#include <vector>

#include "DecaySimulator.h"
#include "Fitter.h"

#include "Minuit2/FCNBase.h"
#include "TGraphErrors.h"
#include "TMatrixD.h"

/*
 * Class representing a simple a + bt + ct^2 polynomial
 */
class SimplePolynomialFunction
{
  public:
    /*
     * Sets the parameters
     */
    SimplePolynomialFunction(double a, double b, double c);

    /*
     * Empty destructor
     */
    ~SimplePolynomialFunction();

    /*
     * Given our parameters, return the function a + bt + ct^2 evaluated at x
     */
    double operator()(double x) const;

  private:
    /*
     * Polynomial coefficients
     */
    double _a{0};
    double _b{0};
    double _c{0};
};

/*
 * Class for calculating chi squared from a SimplePolynomialFunction
 */
class PolynomialFitFcn : public ROOT::Minuit2::FCNBase
{
  public:
    /*
     * Sets the various private members
     */
    PolynomialFitFcn(const std::vector<double>& meas, const std::vector<double>& pos, const std::vector<double>& mvar);

    /*
     * Maybe we need a destructor
     */
    ~PolynomialFitFcn();

    /*
     * Returns _theErrorDef (i still don't really know what this is )
     */
    virtual double Up() const;

    /*
     * Creates a model given a vector of three parameters.
     *
     * Returns the chi squared value between the model given our parameters and the measured data.
     *
     */
    virtual double operator()(const std::vector<double>& parameters) const;

    /*
     * Return the values of the function
     */
    std::vector<double> measurements() const;

    /*
     * Return the positions where our datapoints are evaluated...?
     */
    std::vector<double> positions() const;

    /*
     * Return variances
     */
    std::vector<double> variances() const;

    /*
     * Sets whatever _theErrorDef is
     */
    void setErrorDef(double def);

  private:
    /*
     * Measurements; our datapoints
     */
    std::vector<double> _theMeasurements;

    /*
     * Where our datapoints are evaluated ...?
     */
    std::vector<double> _thePositions;

    /*
     * Errors
     */
    std::vector<double> _theMVariances;

    /*
     * I don't know
     */
    double _theErrorDef;
};

#endif // MINUIT_FITTER_H
