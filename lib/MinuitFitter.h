/*
 * Fit binned data with Minuit
 */
#ifndef MINUIT_FITTER_H
#define MINUIT_FITTER_H

#include <memory>
#include <vector>

#include "DecaySimulator.h"

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
 * Base class for fitting a polynomial using Minuit2
 *
 * Does not contain the necessary operator() method needed for Minuit- child classes should define this.
 */
class BasePolynomialFcn : public ROOT::Minuit2::FCNBase
{
  public:
    /*
     * Sets the various private members
     */
    BasePolynomialFcn(const std::vector<double>& data,
                      const std::vector<double>& times,
                      const std::vector<double>& errors);

    /*
     * Maybe we need a destructor
     */
    ~BasePolynomialFcn();

    /*
     * Returns _theErrorDef (i still don't really know what this is )
     */
    virtual double Up() const;

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

    /*
     * Measurements; our datapoints
     */
    std::vector<double> theMeasurements;

    /*
     * Where our datapoints are evaluated ...?
     */
    std::vector<double> thePositions;

    /*
     * Errors
     */
    std::vector<double> theMVariances;

    /*
     * I don't know
     */
    double theErrorDef;
};

/*
 * Class fitting to a polynomial using chi squared
 */
class PolynomialChiSqFcn : public BasePolynomialFcn
{
  public:
    /*
     * Calls parent constructor
     */
    PolynomialChiSqFcn(const std::vector<double>& data,
                       const std::vector<double>& times,
                       const std::vector<double>& errors);

    ~PolynomialChiSqFcn();

    /*
     * Creates a model given a vector of three parameters.
     *
     * Returns the chi squared value between the model given our parameters and the measured data.
     *
     */
    virtual double operator()(const std::vector<double>& parameters) const;
};

#endif // MINUIT_FITTER_H
