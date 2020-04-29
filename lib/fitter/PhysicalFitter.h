#ifndef PHYSICAL_FITTER_H
#define PHYSICAL_FITTER_H

#include "MinuitScannerBase.h"

/*
 * Class for fitting the ratio of DCS to CF decay times by optimising parameters x, y, rD, Im(Z), Re(Z) and the decay
 * width
 *
 * Expect ratio = rD2 + rD(yRe(Z) + xIm(Z))width*t + (x2 + y2)/4 * (width *t)2
 *
 * Cannot fit without fixing at least one of x, y or a component of Z
 */
class PhysicalFitter : public MinuitScannerBase
{
  public:
    /*
     * Calls parent constructor
     *
     * Must tell this class the options to use when integrating (bin limits etc.), as the model for the fit integrates
     * over the bins to find the expected ratios.
     *
     * Includes option to constrian x and y to their world average values.
     */
    PhysicalFitter(const FitData_t& fitData, const IntegralOptions_t& integralOptions, const bool constrainXY = false);

    /*
     * Fit our data to a second-order polynomial r2 + r(yRZ + xImZ)Gt + (x2+y2)(Gt)2/4 using Minuit2 and the chi-squared
     * method.
     *
     * If xy are not constrained, must fix at least three of {x, y, z_im, z_re, width} before fitting (else fit is
     * poorly defined). If xy are constrained, then only one of {z_im, z_re, width} needs to be fixed.
     *
     * Allocates memory to _plot and _bestFitFunction
     *
     * Populates fitParams
     */
    void fit();

    /*
     * Set this instance's _parameters attribute
     *
     */
    void setPhysicalFitParams(const std::vector<double>& initialParams, const std::vector<double>& initialErrors);

  protected:
    /*
     * Vector of parameters used in the fit
     *
     * Useful for indexing + consistency; we can use this to e.g. fix a parameter by name
     */
    const std::vector<std::string> _paramNames{"x", "y", "r", "z_im", "z_re", "width"};

  private:
    /*
     * Whether to impose the world averages as a constraint onto x and y
     */
    bool _constrainXY = false;
};

#endif // PHYSICAL_FITTER_H
