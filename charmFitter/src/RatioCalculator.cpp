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

/*
RatioCalculator::RatioCalculator(const std::vector<size_t> &denominator,
                                 const std::vector<size_t> &numerator,
                                 const std::vector<double> &binLimits)
{
    // Bin limits should be sorted
    if (!std::is_sorted(binLimits.begin(), binLimits.end())) {
        std::cerr << "Bad time bin limits; should be sorted" << std::endl;
        throw D2K3PiException();
    }
    // binLimits should define every edge of every bin, from the leftmost edge of the smallest bin to the rightmost edge
    // of the highest.
    _binLimits = binLimits;
    _numBins   = binLimits.size() - 1;

    // Find the centres and widths of each time bin + store class attributes
    binCentres = std::vector<double>(_numBins);
    binWidths  = std::vector<double>(_numBins);
    for (size_t i = 0; i < _numBins; ++i) {
        binCentres[i] = 0.5 * (binLimits[i] + binLimits[i + 1]);
        binWidths[i]  = (binLimits[i + 1] - binLimits[i]);
    }

    // Store CF and DCS decay counts as class attributes.
    _denominator = denominator;
    _numerator   = numerator;
}

void RatioCalculator::calculateRatios(void)
{
    // Check that our two datasets are the right length and doesn't contain zeros
    if (std::count(_denominator.begin(), _denominator.end(), 0) ||
        std::count(_numerator.begin(), _numerator.end(), 0)) {
        std::cerr << "Cannot have 0 points in any bin" << std::endl;
        throw D2K3PiException();
    }

    size_t numEntries = _denominator.size();
    if (numEntries != _numerator.size() || numEntries != _numBins) {
        std::cerr
            << "Cannot take ratio of two datasets of different lengths, or with different length to the number of bins."
            << std::endl;
        throw D2K3PiException();
    }

    // Initialise ratio and errors to the right lengths
    ratio = std::vector<double>(_numBins);
    error = std::vector<double>(_numBins);

    // Iterate over the datasets, calculating the ratio of the two datasets and their error
    // Errors in a quotient add in quadrature; assuming our datasets obey Poisson statistics, we find that the
    // fractional error in our ratio is sqrt((N+M)/NM)
    for (size_t i = 0; i < numEntries; ++i) {
        double numeratorVal   = _numerator[i];
        double denominatorVal = _denominator[i];
        double ratioVal       = numeratorVal / denominatorVal;
        ratio[i]              = ratioVal;
        error[i]              = std::sqrt((numeratorVal + denominatorVal) / (numeratorVal * denominatorVal)) * ratioVal;
    }
}
*/

#endif // RATIOCALCULATOR_CPP
