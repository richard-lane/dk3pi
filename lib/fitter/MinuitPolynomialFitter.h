#include "MinuitScannerBase.h"

#ifndef MINUIT_POLY_FITTER_H
#define MINUIT_POLY_FITTER_H

/*
 * Fit to a polynomial (a + bt + ct^2) using the Minuit2 APIs
 */
class MinuitPolynomialFitter : public MinuitScannerBase
{
  public:
    /*
     * Calls parent constructor
     */
    MinuitPolynomialFitter(const FitData_t& fitData);

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
