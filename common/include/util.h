#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <memory>
#include <random>

#include <boost/math/quadrature/trapezoidal.hpp>

#include "D2K3PiError.h"

#include <TCanvas.h>
#include <TLegend.h>
#include <TObject.h>

namespace util
{

/*
 * Parameters that should be passed to TLegend() constructor
 *
 * Co-ordinates are the location of the legend
 * Header is text that is at the top of the legend box (supports ROOTs weird latex thing)
 * Options are the same as for TPave; basically just tells where to draw the box's shadow + whether to draw borders.
 */
typedef struct LegendParams {
    // Positions
    double x1{0};
    double x2{0};
    double y1{0};
    double y2{0};

    // Title
    std::string header = "";

    // Options
    std::string options = "brNDC";
} LegendParams_t;

/*
 * Save a TObject to file
 * The extension in the specified path determines the format of the file.
 */
void saveObjectToFile(TObject *                   myObject,
                      const std::string &         path,
                      const std::string &         drawOptions  = "",
                      const util::LegendParams_t *legendParams = nullptr);

/*
 * Save multiple TObjects to file, useful for e.g. plotting multiple TGraphs on one canvas
 *
 * The objects must be castable to a TObject*; they will be cast back to the type provided
 * in the template parameter before plotting. This is necessary because this function only accepts a vector of TObjects;
 * this is a different type from a vector of e.g. TGraphs!
 *
 * The extension in the specified path determines the format of the file.
 *
 * A legend is mandatory; its parameters are encoded in legendParams
 *
 */
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

/*
 * Plot a vector as a histogram
 *
 * This won't make a publication-quality plot but is useful for checking stuff
 */
void saveHistogram(const std::vector<double> &vector,
                   const std::string &        path,
                   const std::string &        drawOptions = "",
                   const size_t               nBins       = 100);

/*
 * Find which bin a number belongs in
 *
 * 0 for first bin
 */
size_t findBinIndex(const std::vector<double> &binLimits, const double value);

/*
 * Find how many points of myVector belong in each bin defined by binLimits
 *
 */
std::vector<size_t> binVector(const std::vector<double> &myVector, const std::vector<double> &binLimits);

/*
 * Utility function to find bin limits given a dataset and the minimum number of points in each bin
 * Creates as many bins as possible containing exactly the minimum number of points; the last bin may contain more
 * points than the minimum.
 *
 * Note this may result in some very wide bins at the start/end if lowBin/highBin are specified as much lower or higher
 * than the extreme values in the dataset.
 *
 * Params:
 *   dataset         - data to find bins for. Should be sorted.
 *   minPointsPerBin - minimum number of points per bin.
 *   lowBin          - lower edge of the first bin.
 *   highBin         - upper edge of the last bin.
 */
std::vector<double> findBinLimits(const std::vector<double> &dataSet,
                                  const size_t               minPointsPerBin,
                                  const double               lowBin,
                                  const double               highBin);

/*
 * Take two vectors of Re() and Im() values, output a vector of (magnitude, phase) pairs
 * Phase is between 0 and 2pi
 */
std::vector<std::pair<double, double>> reIm2magPhase(const std::vector<double> real,
                                                     const std::vector<double> imaginary);

/*
 * Find mean of container as a double
 */
template <typename Container> double findMean(const Container &c)
{
    double result{0.0};
    for (auto x : c) {
        result += x;
    }
    return result / c.size();
}

/*
 * Find the mean and std dev of a vector
 */
std::pair<double, double> meanAndStdDev(const std::vector<double> &v);

/*
 * Covariance matrix between a vector of datasets
 */
std::vector<std::vector<double>> covarianceMatrix(const std::vector<std::vector<double>> &data);

/*
 * Check if a matrix is square
 */
bool isSquare(const std::vector<std::vector<double>> &matrix);

/*
 * Multiply a vector by a square matrix
 */
std::vector<double> multiply(const std::vector<std::vector<double>> &matrix, const std::vector<double> &vector);

/*
 * Return vectors of random numbers from correlated gaussians
 *
 * Must provide a pointer to a random number generator
 *
 * Returns a vector of vectors of length count
 */
std::vector<std::vector<double>> correlatedGaussianNumbers(const std::shared_ptr<std::mt19937> &   gen,
                                                           const size_t                            count,
                                                           const std::vector<double> &             means,
                                                           const std::vector<std::vector<double>> &covarianceMatrix);

/*
 * Cholesky decomposition of a positive definite symmetric matrix, M = LL^T
 *
 * Doesn't check that the input matrix is positive-definite! No idea what will happen if you pass something else into it
 *
 * Returns lower triangular matrix L
 */
std::vector<std::vector<double>> choleskyDecomp(const std::vector<std::vector<double>> &matrix);

} // namespace util

#endif // UTIL_H
