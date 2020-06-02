/*
 * Utility functions that will be useful in multiple places.
 */
#ifndef UTIL_CPP
#define UTIL_CPP

#include <boost/filesystem.hpp>
#include <iostream>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"

#include "D2K3PiError.h"
#include "util.h"

namespace util
{

void saveObjectToFile(TObject *                   myObject,
                      const std::string &         path,
                      const std::string &         drawOptions,
                      const util::LegendParams_t *legendParams)
{
    TCanvas *c = new TCanvas();
    c->cd();
    myObject->Draw(drawOptions.c_str());

    if (legendParams) {
        TLegend *legend = new TLegend(legendParams->x1,
                                      legendParams->y1,
                                      legendParams->x2,
                                      legendParams->y2,
                                      legendParams->header.c_str(),
                                      legendParams->options.c_str());
        legend->Draw();
        c->SaveAs(path.c_str());
        delete legend;

    } else { // stupid
        c->SaveAs(path.c_str());
    }

    delete c;
}

template <typename T>
void saveObjectsToFile(const std::vector<TObject *> &  myObjects,
                       const std::vector<std::string> &drawOptions,
                       const std::vector<std::string> &legendLabel,
                       const std::string &             path,
                       const LegendParams_t &          legendParams)
{
    // Check that we have the same number of objects as options
    size_t numObjects      = myObjects.size();
    size_t numDrawOptions  = drawOptions.size();
    size_t numLegendLabels = legendLabel.size();
    if (numObjects != numDrawOptions || numLegendLabels != numDrawOptions) {
        std::cerr << "Must pass same number of objects, options and legend labels; have " << numObjects << ", "
                  << numDrawOptions << " and " << numLegendLabels << "." << std::endl;
        throw D2K3PiException();
    }
    if (numObjects == 0) {
        std::cerr << "Cannot plot 0 objects" << std::endl;
        throw D2K3PiException();
    }

    TCanvas *canvas = new TCanvas();
    canvas->cd();

    // Create a legend object; if we're drawing multiple plots on one canvas then we better have a legend to tell them
    // apart
    TLegend *legend = new TLegend(legendParams.x1,
                                  legendParams.y1,
                                  legendParams.x2,
                                  legendParams.y2,
                                  legendParams.header.c_str(),
                                  legendParams.options.c_str());

    for (size_t i = 0; i < numObjects; ++i) {
        T *myObject = (T *)myObjects[i];

        legend->AddEntry(myObject, legendLabel[i].c_str(), "le");
        myObject->Draw(drawOptions[i].c_str());
    }

    legend->Draw();
    canvas->SaveAs(path.c_str());
    delete legend;
    delete canvas;
}

std::vector<size_t> binVector(const std::vector<double> &myVector, const std::vector<double> &binLimits)
{
    size_t              numBins = binLimits.size() - 1;
    std::vector<size_t> numPerBin(numBins);

    // First sort our vector which should make things easier
    std::vector<double> sortedVector = myVector;
    std::sort(sortedVector.begin(), sortedVector.end());

    // Check that our binLimits cover the right range
    double firstPoint = sortedVector[0];
    double lastPoint  = sortedVector.back();
    double lowBin     = binLimits[0];
    double highBin    = binLimits.back();

    if (firstPoint < lowBin || lastPoint > highBin) {
        std::cerr << "bin limits from " << lowBin << " to " << highBin << " do not cover range of data from "
                  << firstPoint << " to " << lastPoint << std::endl;
        throw D2K3PiException();
    }

    // Iterate over the vector
    size_t currentBin = 1;
    for (auto it = sortedVector.begin(); it != sortedVector.end(); ++it) {
        // If this value is more than the current bin limit, we want to put subsequent points in a higher bin.
        while (*it > binLimits[currentBin]) {
            ++currentBin;

            // This error might never get hit but I can't be bothered to think about it right now
            if (currentBin > numBins) {
                std::cerr << "Attempted to insert points into a bin that doesn't exist." << std::endl;
                throw D2K3PiException();
            }
        }

        // currentBin starts at 1 but vector indexing starts at 0.
        ++numPerBin[currentBin - 1];
    }

    return numPerBin;
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
        binLimits[i] = 0.5 * (dataSet[minPointsPerBin * i - 1] + dataSet[minPointsPerBin * i]);
    }

    return binLimits;
}

std::vector<double> exponentialBinLimits(const double maxTime, const double decayConstant, const size_t numBins)
{
    std::vector<double> binLimits{};
    for (size_t i = 0; i <= numBins; ++i) {
        double x = (double)i / numBins;
        double z = 1 - std::exp(-1 * decayConstant * maxTime);
        binLimits.push_back((-1 / decayConstant) * std::log(1 - z * x));
    }

    return binLimits;
}

std::vector<double> expectedParams(const DecayParams_t &phaseSpaceParams)
{
    double expected_a = phaseSpaceParams.r * phaseSpaceParams.r;
    double expected_b = phaseSpaceParams.r *
                        (phaseSpaceParams.y * phaseSpaceParams.z_re + phaseSpaceParams.x * phaseSpaceParams.z_im) *
                        phaseSpaceParams.width;
    double expected_c = 0.25 * (std::pow(phaseSpaceParams.x, 2) + std::pow(phaseSpaceParams.y, 2)) *
                        std::pow(phaseSpaceParams.width, 2);

    return std::vector<double>{expected_a, expected_b, expected_c};
}

std::vector<std::pair<double, double>> reIm2magPhase(const std::vector<double> real,
                                                     const std::vector<double> imaginary)
{
    if (real.size() != imaginary.size()) {
        std::cerr << "im and re vectors different lengths- cannot convert to magnitude + phase" << std::endl;
        throw D2K3PiException();
    }

    std::vector<std::pair<double, double>> outVector(real.size());

    for (size_t i = 0; i < real.size(); ++i) {
        // If our point is (0 + 0i), choose the phase to be 0
        // Otherwise use the normal definition
        if (real[i] == 0. && imaginary[i] == 0.) {
            outVector[i].first  = 0;
            outVector[i].second = 0;
        } else {
            outVector[i].first = std::sqrt(real[i] * real[i] + imaginary[i] * imaginary[i]);
            // Use atan2 to determine the quadrant
            double phase = std::fmod(std::atan2(imaginary[i], real[i]), 2 * M_PI);
            if (phase < 0) {
                phase += 2 * M_PI;
            }
            outVector[i].second = phase;
        }
    }
    return outVector;
}

// Explicitly instantiate the types we want to use for saving objects to file
// This is Ugly and Bad but also i dont mind
template void saveObjectsToFile<TGraph>(const std::vector<TObject *> &  myObjects,
                                        const std::vector<std::string> &drawOptions,
                                        const std::vector<std::string> &legendLabel,
                                        const std::string &             path,
                                        const LegendParams_t &          legendParams);

template void saveObjectsToFile<TGraphErrors>(const std::vector<TObject *> &  myObjects,
                                              const std::vector<std::string> &drawOptions,
                                              const std::vector<std::string> &legendLabel,
                                              const std::string &             path,
                                              const LegendParams_t &          legendParams);

} // namespace util

#endif // UTIL_CPP
