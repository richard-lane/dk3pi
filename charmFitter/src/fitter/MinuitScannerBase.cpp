
#include <iostream>

#include <boost/progress.hpp>

#include "D2K3PiError.h"
#include "MinuitScannerBase.h"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnScan.h"

MinuitScannerBase::MinuitScannerBase(const FitData_t& fitData)
    : MinuitFitterBase(fitData)
{
    ;
}

std::vector<std::pair<double, double>>
MinuitScannerBase::chiSqParameterScan(const size_t i, const size_t numPoints, const double low, const double high)
{
    // Find a vector of the parameter values we're interested in
    // This already gets done in the _scanParameter function but i dont care
    std::vector<double> parameterVals(numPoints, 0.0);
    double              step = (high - low) / (numPoints - 1);
    for (size_t k = 0; k < numPoints; ++k) {
        parameterVals[k] = low + k * step;
    }

    // Scan the parameter and populate parameterScan
    std::vector<std::pair<double, double>> parameterScan(numPoints);
    std::vector<double>                    chiSqValues = _scanParameter(i, numPoints, low, high);
    for (size_t k = 0; k < numPoints; ++k) {
        parameterScan[k] = std::make_pair(parameterVals[k], chiSqValues[k]);
    }

    return parameterScan;
}

std::vector<std::vector<double>> MinuitScannerBase::twoDParamScan(const size_t                    i,
                                                                  const size_t                    j,
                                                                  const size_t                    iPoints,
                                                                  const size_t                    jPoints,
                                                                  const std::pair<double, double> iRange,
                                                                  const std::pair<double, double> jRange)
{
    if (i >= _parameters->Params().size() || j >= _parameters->Params().size()) {
        std::cerr << "Cannot scan params " << i << ", " << j << "; only have " << fitParams.fitParams.size()
                  << " params." << std::endl;
        throw D2K3PiException();
    }

    if (i == j) {
        std::cerr << "Cannot perfom 2d scan on a parameter against itself" << std::endl;
        throw D2K3PiException();
    }

    // Create vectors to store the i, j values we're interested in, and a vector to store our result in
    std::vector<double>              iVals(iPoints);
    std::vector<double>              jVals(jPoints);
    std::vector<std::vector<double>> twoDParameterScan(iPoints * jPoints);

    // Fill vectors of i and j values
    // Technically we could do all this in the big loop where we do the minimisation but i think this is clearer to
    // think about
    double iStep = (iRange.second - iRange.first) / (iPoints - 1);
    double jStep = (jRange.second - jRange.first) / (jPoints - 1);
    for (size_t k = 0; k < iPoints; ++k) {
        iVals[k] = iRange.first + k * iStep;
    }
    for (size_t k = 0; k < jPoints; ++k) {
        jVals[k] = jRange.first + k * jStep;
    }

    // Loop over our j Values, performing a scan over i  for each one
    _parameters->Fix(j);
    boost::progress_display showProgress(jPoints);
    for (size_t jIndex = 0; jIndex < jPoints; ++jIndex) {
        _parameters->SetValue(j, jVals[jIndex]);
        std::vector<double> values = _scanParameter(i, iPoints, iRange.first, iRange.second);
        for (size_t n = 0; n < values.size(); ++n) {
            twoDParameterScan[jIndex + jPoints * n] = std::vector<double>{iVals[n], jVals[jIndex], values[n]};
        }
        ++showProgress;
    }
    _parameters->Release(j);

    return twoDParameterScan;
}

std::vector<double>
MinuitScannerBase::_scanParameter(const size_t i, const size_t numPoints, const double low, const double high)
{
    // Check that we have sensible low and high
    if (low > high) {
        std::cerr << "Low value must be below high value for parameter scan." << std::endl;
        throw D2K3PiException();
    }

    // Check that our parameter index is in range
    if (i >= _parameters->Params().size()) {
        std::cerr << "Cannot scan parameter " << i << "; out of range" << std::endl;
        throw D2K3PiException();
    }

    // Create a vector of the right length to store our function values
    std::vector<double> values(numPoints);

    // Find a vector of the parameter values we're interested in
    std::vector<double> parameterVals(numPoints, 0.0);
    double              step = (high - low) / (numPoints - 1);
    for (size_t k = 0; k < numPoints; ++k) {
        parameterVals[k] = low + k * step;
    }

    // For every value of the parameter we're interested in, fix the param, perform a fit + populate values
    _parameters->Fix(i);
    for (size_t k = 0; k < numPoints; ++k) {
        _parameters->SetValue(i, parameterVals[k]);

        // Perform a fit and insert the statistic value we get into values[]
        fit();
        values[k] = fitParams.fitStatistic;
    }
    _parameters->Release(i);

    return values;
}