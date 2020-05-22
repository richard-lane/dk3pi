#ifndef RATIOCALCULATOR_CPP
#define RATIOCALCULATOR_CPP

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "D2K3PiError.h"
#include "RatioCalculator.h"
#include "util.h"

#include "TMath.h"

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

#endif // RATIOCALCULATOR_CPP
