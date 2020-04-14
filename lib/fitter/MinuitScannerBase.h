#ifndef MINUIT_SCANNER_BASE_H
#define MINUIT_SCANNER_BASE_H

#include <utility>
#include <vector>

#include "MinuitFitterBase.h"

class MinuitScannerBase : public MinuitFitterBase
{
  public:
    /*
     * Scan the ith parameter as defined in fitParams.fitParams
     *
     * Cannot have more than 100 points due to some of Minuit's internal limitations.
     * By default runs a scan from +-2sigma, but can optionally be specified by setting low and high.
     *
     * Populates parameterScan
     */
    void chiSqParameterScan(const size_t i, const size_t numPoints, const double low, const double high);

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

  protected:
    /*
     * Calls parent constructor
     */
    MinuitScannerBase(const FitData_t& fitData);

  private:
    /*
     * Perform a series of fits, holding parameter i fixed between the values specified in low and high.
     *
     * Returns a vector of the minimised Minuit2 function values
     */
    std::vector<double> _scanParameter(const size_t i, const size_t numPoints, const double low, const double high);
};

#endif // MINUIT_SCANNER_BASE_H
