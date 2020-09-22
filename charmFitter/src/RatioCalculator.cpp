#ifndef RATIOCALCULATOR_CPP
#define RATIOCALCULATOR_CPP

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "RatioCalculator.h"
#include "util.h"

#include "TMath.h"

namespace RatioCalculator
{

std::pair<std::vector<double>, std::vector<double>> ratioAndError(const std::vector<size_t> &denominator,
                                                                  const std::vector<size_t> &numerator)
{
    // Check that our two datasets are the right length and doesn't contain zeros
    if (std::count(denominator.begin(), denominator.end(), 0) || std::count(numerator.begin(), numerator.end(), 0)) {
        throw ZerosInData();
    }
    size_t numBins = denominator.size();
    if (numBins != numerator.size()) {
        throw WrongSizeDataset();
    }

    // Initialise ratio and errors to the right lengths
    std::vector<double> ratio(numBins);
    std::vector<double> error(numBins);

    // Iterate over the datasets, calculating the ratio of the two datasets and their error
    // Errors in a quotient add in quadrature; assuming our datasets obey Poisson statistics, we find that the
    // fractional error in our ratio is sqrt((N+M)/NM)
    for (size_t i = 0; i < numBins; ++i) {
        ratio[i] = static_cast<double>(numerator[i]) / denominator[i];
        error[i] = std::sqrt((static_cast<double>(numerator[i]) + denominator[i]) / (numerator[i] * denominator[i])) *
                   ratio[i];
    }

    return std::pair<std::vector<double>, std::vector<double>>(ratio, error);
}

} // namespace RatioCalculator

#endif // RATIOCALCULATOR_CPP
