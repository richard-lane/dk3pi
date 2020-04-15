#ifndef MINUIT_SCANNER_BASE_H
#define MINUIT_SCANNER_BASE_H

#include <utility>
#include <vector>

#include "MinuitFitterBase.h"

class MinuitScannerBase : public MinuitFitterBase
{
  public:
    /*
     * Scan the ith parameter as defined in _parameters
     *
     * Also may pass in a default value for chi squared; if this is set, then this will be used if a fit fails (a
     * warning message will be emitted).
     *
     * Populates parameterScan
     */
    void chiSqParameterScan(const size_t  i,
                            const size_t  numPoints,
                            const double  low,
                            const double  high,
                            const double* defaultChiSq = nullptr);

    /*
     * Scan the i and jth parameters between the specified limits
     *
     * Also may pass in a default value for chi squared; if this is set, then this will be used if a fit fails (a
     * warning message will be emitted).
     *
     * Populates twoDParameterScan
     *
     *
     * TODO this has too many args
     */
    void twoDParamScan(const size_t  i,
                       const size_t  j,
                       const size_t  iPoints,
                       const size_t  jPoints,
                       const double  iLow,
                       const double  iHigh,
                       const double  jLow,
                       const double  jHigh,
                       const double* defaultChiSq = nullptr);

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
     * Also fixes the parameters specified in additionalFixParams (can be empty)
     *
     * If a default value of chi squared is provided, failed fits will be ignored + the value given will be used instead
     * (a warning message is still emitted)
     *
     * Returns a vector of the minimised Minuit2 function values
     */
    std::vector<double> _scanParameter(const size_t               i,
                                       const size_t               numPoints,
                                       const double               low,
                                       const double               high,
                                       const std::vector<size_t>& additionalFixParams,
                                       const double*              defaultChiSq = nullptr);
};

#endif // MINUIT_SCANNER_BASE_H
