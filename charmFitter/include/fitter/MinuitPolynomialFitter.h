#ifndef MINUIT_POLY_FITTER_H
#define MINUIT_POLY_FITTER_H

#include "MinuitScannerBase.h"

/*
 * Class fitting to a polynomial using chi squared
 *
 * Fits by integrating over the bins
 */
class PolynomialChiSqFcn : public MyBaseFcn
{
  public:
    /*
     * Calls parent constructor
     */
    PolynomialChiSqFcn(const std::vector<double>& data,
                       const std::vector<double>& times,
                       const std::vector<double>& errors,
                       const IntegralOptions_t&   integralOptions);

    ~PolynomialChiSqFcn();

    /*
     * Creates a model given a vector of three parameters.
     *
     * Returns the chi squared value between the model given our parameters and the measured data.
     *
     */
    virtual double operator()(const std::vector<double>& parameters) const;

  private:
    /*
     * options to use while integrating
     */
    IntegralOptions_t _integralOptions;
};

/*
 * Class fitting to a polynomial using chi squared
 *
 * Does not integrate over the bins
 */
class PolynomialChiSqFcnNoIntegral : public MyBaseFcn
{
  public:
    /*
     * Calls parent constructor
     */
    PolynomialChiSqFcnNoIntegral(const std::vector<double>& data,
                                 const std::vector<double>& times,
                                 const std::vector<double>& errors);

    ~PolynomialChiSqFcnNoIntegral();

    /*
     * Creates a model given a vector of three parameters.
     *
     * Returns the chi squared value between the model given our parameters and the measured data.
     *
     */
    virtual double operator()(const std::vector<double>& parameters) const;
};

/*
 * Fit to a polynomial (a + bt + ct^2) using the Minuit2 APIs
 */
class MinuitPolynomialFitter : public MinuitScannerBase
{
  public:
    /*
     * Calls parent constructor
     *
     * If integralOptions.integrate == true , the fitter will integrate over each bin when calculating chi squared. This
     * will make the fit way slower but also maybe more accurate?
     */
    MinuitPolynomialFitter(const FitData_t& fitData, const IntegralOptions_t& integralOptions);

    /*
     * Fit our data to a second-order polynomial a + bt + ct^2 using Minuit2 and the chi-squared method.
     *
     * User should call setPolynomialParams first to set an initial guess at the parameters and their errors
     *
     * Allocates memory to _plot and _bestFitFunction
     *
     * Populates fitParams
     */
    void fit();

    /*
     * Set initial values + errors for a, b, and c
     */
    void setPolynomialParams(const std::vector<double>& initialParams, const std::vector<double>& initialErrors);
};

#endif // MINUIT_POLY_FITTER_H
