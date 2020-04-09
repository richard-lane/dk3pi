/*
 * Fit binned data
 */
#ifndef FITTER_H
#define FITTER_H

#include <memory>
#include <vector>

#include "BaseFitter.h"
#include "DecaySimulator.h"
#include "FitterUtils.h"
#include "MinuitFitter.h"
#include "MinuitFitterBase.h"
#include "MinuitPolynomialFitter.h"
#include "util.h"

#include "Minuit2/FunctionMinimum.h"
#include "TGraphErrors.h"
#include "TMatrixD.h"

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

/*
 * Class for perfoming scans of the physical fit parameters x, y etc.
 */
class ParamScanner : public PhysicalFitter
{
  public:
    /*
     * Calls parent constructor
     */
    ParamScanner(const FitData_t& fitData);

    /*
     * Scan the ith parameter as defined in fitParams.fitParams
     *
     * Cannot have more than 100 points due to some of Minuit's internal limitations.
     * By default runs a scan from +-2sigma, but can optionally be specified by setting low and high.
     *
     * Populates parameterScan
     */
    void chiSqParameterScan(const size_t i, const size_t numPoints, const double low = 0., const double high = 0.);

    /*
     * Scan the i and jth parameters between the specified limits
     *
     * Populates twoDParameterScan
     */
    void twoDParamScan(const size_t i,
                       const size_t j,
                       const size_t iPoints,
                       const size_t jPoints,
                       const double iLow,
                       const double iHigh,
                       const double jLow,
                       const double jHigh);

    /*
     * Vector of pairs describing a parameter scan
     */
    std::vector<std::pair<double, double>> parameterScan;

    /*
     * Vector of tuples describing a 2d parameter scan
     *
     * Scans parameters i and j to find chi squared values; result is a vector of (i_value, j_value, chi_squared)
     */
    std::vector<std::vector<double>> twoDParameterScan;
};

#endif // FITTER_H
