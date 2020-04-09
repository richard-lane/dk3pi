#ifndef MINUIT_POLY_SCAN_H
#define MINUIT_POLY_SCAN_H

#include "MinuitPolynomialFitter.h"

/*
 * Class for perfoming scans of the polynomial fit parameters a, b, c
 */
class MinuitPolyScan : public MinuitPolynomialFitter
{
  public:
    /*
     * Calls parent constructor
     */
    MinuitPolyScan(const FitData_t& fitData);

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

#endif // MINUIT_POLY_SCAN_H