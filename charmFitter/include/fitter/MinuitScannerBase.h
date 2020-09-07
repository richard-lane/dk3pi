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
     * Returns a vector std::pairs (parameter val, fit statistic val)
     */
    std::vector<std::pair<double, double>>
    chiSqParameterScan(const size_t i, const size_t numPoints, const double low, const double high);

    /*
     * Scan the i and jth parameters between the specified limits (iRange) and (jRange)
     *
     * Populates twoDParameterScan
     *
     * returns a vector of (parameter i val, parameter j val, fit statistic value)
     *
     */
    std::vector<std::vector<double>> twoDParamScan(const size_t                    i,
                                                   const size_t                    j,
                                                   const size_t                    iPoints,
                                                   const size_t                    jPoints,
                                                   const std::pair<double, double> iRange,
                                                   const std::pair<double, double> jRange);

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
