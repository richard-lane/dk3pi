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
#include "PhysicalFitter.h"
#include "util.h"

#include "Minuit2/FunctionMinimum.h"
#include "TGraphErrors.h"
#include "TMatrixD.h"

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
