#include "MinuitFitterBase.h"

#ifndef MINUIT_POLY_FITTER_H
#define MINUIT_POLY_FITTER_H

/*
 * Fit to a polynomial (a + bt + ct^2) using the Minuit2 APIs
 */
class MinuitPolynomialFitter : public MinuitFitterBase
{
  public:
    /*
     * Calls parent constructor
     */
    MinuitPolynomialFitter(const FitData_t& fitData);

    /*
     * Fit our data to a second-order polynomial a + bt + ct^2 using Minuit2 and the chi-squared method.
     *
     * The user should provide an initial guess at the parameters and their errors
     * Parameters are {x, y, r, z_im, z_re, width}
     *
     * FitMethod tells the fitter whether to use chi squared or maximum likelihood (max likelihood isn't actually
     * implemented)
     *
     * Allocates memory to _plot and _bestFitPlot
     *
     * Populates fitParams
     */
    void fit(const std::vector<double>& initialParams,
             const std::vector<double>& initialErrors,
             const FitAlgorithm_t&      FitMethod);
};

#endif // MINUIT_POLY_FITTER_H
