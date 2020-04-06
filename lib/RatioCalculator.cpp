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

RatioCalculator::RatioCalculator(const std::vector<double> &cfDecayTimes,
                                 const std::vector<double> &dcsDecayTimes,
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

    // Find the centres and widths of each time bin
    binCentres.assign(_numBins, -1);
    binWidths.assign(_numBins, -1);
    binCentres.shrink_to_fit();
    binCentres.shrink_to_fit();

    for (size_t i = 0; i < _numBins; ++i) {
        binCentres[i] = 0.5 * (binLimits[i] + binLimits[i + 1]);
        binWidths[i]  = (binLimits[i + 1] - binLimits[i]);
    }

    // Store CF and DCS decay times as class attributes.
    _cfDecayTimes  = cfDecayTimes;
    _dcsDecayTimes = dcsDecayTimes;
}

std::vector<std::pair<double, double>> RatioCalculator::findRatioAndError(const std::vector<size_t> &numerator,
                                                                          const std::vector<size_t> &denominator)
{
    std::vector<std::pair<double, double>> ratiosAndErrors(_numBins);

    // Check that our two datasets are the right length
    size_t numEntries = denominator.size();
    if (numEntries != numerator.size() || numEntries != _numBins) {
        std::cerr
            << "Cannot take ratio of two datasets of different lengths, or with different length to the number of bins."
            << std::endl;
        throw D2K3PiException();
    }

    // Iterate over the datasets, calculating the ratio of the two datasets and their error
    // Errors in a quotient add in quadrature; assuming our datasets obey Poisson statistics, we find that the
    // fractional error in our ratio is sqrt((N+M)/NM)
    for (size_t i = 0; i < numEntries; ++i) {
        double numeratorVal      = numerator[i];
        double denominatorVal    = denominator[i];
        double ratio             = numeratorVal / denominatorVal;
        ratiosAndErrors[i].first = ratio;
        ratiosAndErrors[i].second =
            std::sqrt((numeratorVal + denominatorVal) / (numeratorVal * denominatorVal)) * ratio;
    }

    return ratiosAndErrors;
}

void RatioCalculator::findNumPointsPerBin(const std::string &path)
{
    std::ofstream f(path.c_str());
    f << "CF\tDCS" << std::endl;

    for (size_t i = 0; i < _numBins; ++i) {
        f << std::to_string(_numCfPerBin[i]) << "\t" << std::to_string(_numDcsPerBin[i]) << std::endl;
    }
}

void RatioCalculator::calculateRatios(void)
{
    // Bin both the CF and DCS vectors
    _numCfPerBin  = util::binVector(_cfDecayTimes, _binLimits);
    _numDcsPerBin = util::binVector(_dcsDecayTimes, _binLimits);

    if (std::count(_numDcsPerBin.begin(), _numDcsPerBin.end(), 0) ||
        std::count(_numCfPerBin.begin(), _numCfPerBin.end(), 0)) {
        std::cerr << "Cannot have 0 points in any bin" << std::endl;
        throw D2K3PiException();
    }

    // Find the ratio of CF to DCS vectors and populate the relevant attributes.
    std::vector<std::pair<double, double>> ratiosAndErrors = findRatioAndError(_numDcsPerBin, _numCfPerBin);

    // Assign the ratios and errors to the right attributes.
    ratio.resize(_numBins, 0);
    error.resize(_numBins, 0);

    // Set ratios and errors
    for (size_t i = 0; i < _numBins; ++i) {
        ratio[i] = ratiosAndErrors[i].first;
        error[i] = ratiosAndErrors[i].second;
    }
}

#endif // RATIOCALCULATOR_CPP
