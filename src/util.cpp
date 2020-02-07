/*
 * Utility functions that will be useful in multiple places.
 */
#ifndef UTIL_CPP
#define UTIL_CPP

#include <boost/filesystem.hpp>
#include <iostream>

#include "TCanvas.h"

#include "../include/D2K3PiError.h"
#include "../include/util.h"

namespace util
{

boost::filesystem::path concatPaths(std::string plotDir, std::string plotName, std::string fileExtension)
{
    boost::filesystem::path dir(plotDir);
    boost::filesystem::path file(plotName + fileExtension);
    return dir / file;
}

void saveToFile(TObject *myObject, const std::string &path, const std::string &drawOptions)
{
    TCanvas *c = new TCanvas();
    c->cd();
    myObject->Draw(drawOptions.c_str());
    c->SaveAs(path.c_str());

    delete c;
}

std::vector<double> findBinLimits(const std::vector<double> &dataSet,
                                  const size_t               minPointsPerBin,
                                  const double               lowBin,
                                  const double               highBin)
{
    // Check that our dataset is sorted
    if (!std::is_sorted(dataSet.begin(), dataSet.end())) {
        std::cerr << "Data set must be sorted before bin limits are found" << std::endl;
        throw D2K3PiException();
    }

    // Check we haven't been passed 0 for minPointsPerBin
    if (minPointsPerBin == 0) {
        std::cerr << "Cannot have a minimum of 0 points per bin" << std::endl;
        throw D2K3PiException();
    }

    // Check our dataset has more points than the minimum points per bin
    size_t numPoints = dataSet.size();
    if (numPoints < minPointsPerBin) {
        std::cerr << "Cannot sort a dataset of length " << numPoints << " into bins of size " << minPointsPerBin
                  << std::endl;
        throw D2K3PiException();
    }

    // Check our lower and higher bin limits cover the entire dataset
    if (lowBin >= dataSet[0] || highBin <= dataSet[numPoints - 1]) {
        std::cerr << "Range of points (" << lowBin << ", " << highBin << ") does not cover entire dataset."
                  << std::endl;
        throw D2K3PiException();
    }

    // Find how many bins to use. We want to use as few bins as possible so that we have ~minPointsPerBin in each bin
    size_t numBins = numPoints / minPointsPerBin;

    // Find where the bin edges should be placed
    std::vector<double> binLimits(numBins + 1, NAN);
    binLimits[0]       = lowBin;
    binLimits[numBins] = highBin;
    for (size_t i = 1; i < numBins; ++i) {
        // Place our interim bin limits halfway between a point and its neighbour
        binLimits[i] = 0.5 * (dataSet[numBins * i - 1] + (dataSet[numBins * i]));
    }

    return binLimits;
}

} // namespace util

#endif // UTIL_CPP
