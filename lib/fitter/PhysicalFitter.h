#ifndef PHYSICAL_FITTER_H
#define PHYSICAL_FITTER_H

#include "MinuitFitterBase.h"

/*
 * Class for fitting the ratio of DCS to CF decay times by optimising parameters x, y, rD, Im(Z), Re(Z) and the decay
 * width
 *
 * Expect ratio = rD2 + rD(yRe(Z) + xIm(Z))width*t + (x2 + y2)/4 * (width *t)2
 *
 * Cannot fit without fixing at least one of x, y or a component of Z
 */
class PhysicalFitter : public MinuitFitter
{
  public:
    /*
     * Calls parent constructor
     */
    PhysicalFitter(const FitData_t& fitData);

    /*
     * Fit our data to a second-order polynomial r2 + r(yRZ + xImZ)Gt + (x2+y2)(Gt)2/4 using Minuit2 and the chi-squared
     * method.
     *
     * The user should provide an initial guess at the parameters and their errors
     * Also provide a vector of parameters to fix; parameter numbering is defined by _paramNames (should not be empty).
     * e.g. fix parameters 2 and 3 to 0.2, 0.3 respectively by passing {(2, 0.2), (3, 0,3)}
     *
     * Allocates memory to _plot and _bestFitPlot
     *
     * Populates fitParams
     */
    void fit(const std::vector<double>&                    initialParams,
             const std::vector<double>&                    initialErrors,
             const FitAlgorithm_t&                         FitMethod,
             const std::vector<std::pair<size_t, double>>& fixParams);

  protected:
    /*
     * Vector of parameters used in the fit
     *
     * Useful for indexing + consistency; we can use this to e.g. fix a parameter by name
     */
    const std::vector<std::string> _paramNames{"x", "y", "r", "z_im", "z_re", "width"};
};

#endif // PHYSICAL_FITTER_H
