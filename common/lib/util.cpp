/*
 * Utility functions that will be useful in multiple places.
 */
#include <boost/filesystem.hpp>
#include <cmath>
#include <iostream>

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>

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
    canvas->SetLeftMargin(0.15); // Magic number means you can actually read axis labels

    // Create a legend object; if we're drawing multiple plots on one canvas then we better have a legend to tell them
    // apart
    TLegend *legend = new TLegend(legendParams.x1,
                                  legendParams.y1,
                                  legendParams.x2,
                                  legendParams.y2,
                                  legendParams.header.c_str(),
                                  legendParams.options.c_str());
    legend->SetTextSize(0.03);

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

void saveHistogram(const std::vector<double> &vector,
                   const std::string &        path,
                   const std::string &        drawOptions,
                   const size_t               nBins)
{
    double min = 0.9 * *std::min_element(vector.begin(), vector.end());
    double max = 1.1 * *std::max_element(vector.begin(), vector.end());

    std::unique_ptr<TH1> hist(new TH1D(path.c_str(), path.c_str(), nBins, min, max));
    hist->FillN(vector.size(), vector.data(), nullptr);
    saveObjectToFile(hist.get(), path, drawOptions);
}

size_t findBinIndex(const std::vector<double> &binLimits, const double value)
{
    // Check that the value fits in the bin range
    double lowBin  = binLimits[0];
    double highBin = binLimits.back();
    if (value < lowBin || value > highBin) {
        std::cerr << "bin limits from " << lowBin << " to " << highBin << " do not cover range of data from " << lowBin
                  << " to " << highBin << std::endl;
        throw D2K3PiException();
    }
    size_t currentBin = 1;
    while (value > binLimits[currentBin]) {
        ++currentBin;
    }

    return currentBin - 1;
}

std::vector<size_t> binVector(const std::vector<double> &myVector, const std::vector<double> &binLimits)
{
    size_t              numBins = binLimits.size() - 1;
    std::vector<size_t> numPerBin(numBins);

    // Iterate over the vector, placing each point into a bin and incrementing that bin's value
    // This did something cleverer before where I sorted the vector and binned the sorted vector in order, but maybe
    // this is nicer
    for (auto it = myVector.begin(); it != myVector.end(); ++it) {
        ++numPerBin[findBinIndex(binLimits, *it)];
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

std::pair<double, double> meanAndStdDev(const std::vector<double> &v)
{
    double mean = findMean(v);

    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(), std::bind2nd(std::minus<double>(), mean));
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev  = std::sqrt(sq_sum / v.size());

    return std::make_pair(mean, stdev);
}

std::vector<std::vector<double>> covarianceMatrix(const std::vector<std::vector<double>> &data)
{
    size_t numDatasets = data.size();
    size_t dataLength  = data[0].size();

    // Check all datasets are the same length
    for (auto it = data.begin() + 1; it != data.end(); ++it) {
        if (it->size() != dataLength) {
            std::cerr << "passed in datasets of different lengths; cannot calculate covariance" << std::endl;
            throw D2K3PiException();
        }
    }

    // Initialise cov matrix
    std::vector<std::vector<double>> outMatrix(
        numDatasets, std::vector<double>(numDatasets, std::numeric_limits<double>::quiet_NaN()));

    // Fill it
    for (size_t i = 0; i < numDatasets; ++i) {
        double iMean = meanAndStdDev(data[i]).first;

        for (size_t j = 0; j < numDatasets; ++j) {
            double jMean = meanAndStdDev(data[j]).first;
            double cov{0.0};
            for (size_t k = 0; k < dataLength; ++k) {
                cov += (data[i][k] - iMean) * (data[j][k] - jMean);
            }
            outMatrix[i][j] = cov / dataLength;
        }
    }
    return outMatrix;
}

bool isSquare(const std::vector<std::vector<double>> &matrix)
{
    size_t matrixDimension = matrix.size();
    for (size_t i = 0; i < matrixDimension; ++i) {
        if (matrix[i].size() != matrixDimension) {
            return false;
        }
    }
    return true;
}

std::vector<double> multiply(const std::vector<std::vector<double>> &matrix, const std::vector<double> &vector)
{
    if (!isSquare(matrix)) {
        std::cerr << "matrix not square; cannot perform Cholesky decomposition" << std::endl;
        throw D2K3PiException();
    }

    if (matrix.size() != vector.size()) {
        std::cerr << "matrix and vector incompatible sizes" << std::endl;
        throw D2K3PiException();
    }

    std::vector<double> outVector(vector.size(), 0);

    for (size_t i = 0; i < vector.size(); ++i) {
        for (size_t j = 0; j < vector.size(); ++j) {
            outVector[i] += matrix[i][j] * vector[j];
        }
    }

    return outVector;
}

std::vector<std::vector<double>> correlatedGaussianNumbers(const std::shared_ptr<std::mt19937> &   gen,
                                                           const size_t                            count,
                                                           const std::vector<double> &             means,
                                                           const std::vector<std::vector<double>> &covarianceMatrix)
{
    size_t dimension = covarianceMatrix.size();
    if (means.size() != dimension) {
        std::cerr << "incompatible means + covariance matrix provided" << std::endl;
        throw D2K3PiException();
    }

    // Gaussian generator
    std::normal_distribution<double> distribution(0, 1);

    // Find Cholesky decomposition of our covariance matrix
    std::vector<std::vector<double>> lowerTriangular = choleskyDecomp(covarianceMatrix);

    // Create vectors of normally-distributed random numbers
    std::vector<std::vector<double>> randomNumbers(dimension, std::vector<double>(count, -1));
    for (size_t i = 0; i < dimension; ++i) {
        for (size_t j = 0; j < count; ++j) {
            randomNumbers[i][j] = distribution(*gen);
        }
    }

    // Transform them to correlated random numbers
    for (size_t i = 0; i < count; ++i) {
        // Create a vector of uncorrelated random numbers
        std::vector<double> uncorrelated(dimension, -1);
        for (size_t j = 0; j < dimension; ++j) {
            uncorrelated[j] = randomNumbers[j][i];
        }

        // multiply this vector by the lower triangular matrix to get correlated random numbers
        uncorrelated = multiply(lowerTriangular, uncorrelated);

        // Insert them back into the vector of vectors, giving them the correct mean
        for (size_t j = 0; j < dimension; ++j) {
            randomNumbers[j][i] = means[j] + uncorrelated[j];
        }
    }

    return randomNumbers;
}

std::vector<std::vector<double>> choleskyDecomp(const std::vector<std::vector<double>> &matrix)
{
    size_t matrixDimension = matrix.size();
    if (!isSquare(matrix)) {
        std::cerr << "matrix not square; cannot perform Cholesky decomposition" << std::endl;
        throw D2K3PiException();
    }

    // Check that it is symmetric
    for (size_t i = 0; i < matrixDimension; ++i) {
        for (size_t j = 0; j < matrixDimension; ++j) {
            // No the most efficient check but should be fine
            if (matrix[i][j] != matrix[j][i]) {
                std::cerr << "matrix not symmetric; cannot perform Cholesky decomposition" << std::endl;
                throw D2K3PiException();
            }
        }
    }

    // Initialise our output matrix to NaN
    std::vector<std::vector<double>> lower(
        matrixDimension, std::vector<double>(matrixDimension, std::numeric_limits<double>::quiet_NaN()));

    // Cholesky-Branachiewicz algorithm
    for (size_t i = 0; i < matrixDimension; ++i) {
        for (size_t j = 0; j < matrixDimension; ++j) {

            // Upper triangle
            if (j > i) {
                lower[i][j] = 0;
            }

            // Diagonal
            else if (i == j) {
                double sum{0.0};
                for (int k = 0; k < (int)j; ++k) {
                    sum += lower[j][k] * lower[j][k];
                }
                lower[j][j] = std::sqrt(matrix[j][j] - sum);
            }

            // Lower triangle
            else {
                double sum{0.0};
                for (int k = 0; k < int(j); ++k) {
                    sum += lower[i][k] * lower[j][k];
                }
                lower[i][j] = (matrix[i][j] - sum) / lower[j][j];
            }
        }
    }
    return lower;
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
