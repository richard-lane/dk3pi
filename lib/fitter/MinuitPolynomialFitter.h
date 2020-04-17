#ifndef MINUIT_POLY_FITTER_H
#define MINUIT_POLY_FITTER_H

#include "MinuitFitter.h"
#include "MinuitScannerBase.h"

/*
 * Fit to a polynomial (a + bt + ct^2) using the Minuit2 APIs
 */
class MinuitPolynomialFitter : public MinuitScannerBase
{
  public:
    /*
     * Calls parent constructor
     *
     * This class needs to know about bin limits and decay width, as the fit is performed by integrating the expected
     * functions over each bin
     *
     * If integralOptions are provided, the fitter will integrate over each bin when calculating chi squared. This will
     * make the fit way slower, but also maybe more accurate?
     */
    MinuitPolynomialFitter(const FitData_t&           fitData,
                           const std::vector<double>& binLimits,
                           const double               width,
                           const IntegralOptions_t*   integralOptions = nullptr);

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
