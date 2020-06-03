#ifndef UTIL_H
#define UTIL_H

#include <random>

#include <boost/filesystem.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

#include "D2K3PiError.h"

#include <TObject.h>

/*
 * Struct encapsulating the parameters needed to simulate decays.
 * Mixing params x, y
 * Phase space params r, Im(Z) and Re(Z)
 * Particle data \Gamma for the D meson
 *
 * Parameters are all initalised to zero by default and should be set by the user.
 *
 * Should really be in the util namespace
 */
typedef struct DecayParameters {
    // Mixing params
    double x{0.0};
    double y{0.0};

    // Phase space params
    double r{0.0};
    double z_im{0.0};
    double z_re{0.0};

    // Particle data
    double width{0.0};
} DecayParams_t;

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
                       const LegendParams_t &          legendParams);

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
 * Find exponentially-distributed limits of numBins bins from 0 up to maxTime, whose width is described by
 * exp(-decayConstant * time)
 */
std::vector<double> exponentialBinLimits(const double maxTime, const double decayConstant, const size_t numBins);

/*
 * Find our expected a, b and c in (a + bt + ct^2) from a set of phase space parameters.
 *
 * Returns a vector of {a, b, c}
 */
std::vector<double> expectedParams(const DecayParams_t &phaseSpaceParams);

/*
 * Take two vectors of Re() and Im() values, output a vector of (magnitude, phase) pairs
 * Phase is between 0 and 2pi
 */
std::vector<std::pair<double, double>> reIm2magPhase(const std::vector<double> real,
                                                     const std::vector<double> imaginary);

/*
 * Use adaptive quadrature to find an approximation to the integral of f between low and high limits
 */
template <typename Func>
double adadptiveTrapQuad(Func         f,
                         const double low,
                         const double high,
                         const double tolerance,
                         const size_t maxRefinements,
                         double *     errorEstimate = nullptr)
{
    return boost::math::quadrature::trapezoidal(f, low, high, tolerance, maxRefinements, errorEstimate);
}

/*
 * Use Gauss-Legendre quadrature to find an approximation to the integral of f between low and high limits
 *
 * Evaluates the function at 15 points, as the weights and abcissa have been precalculated for this number of points.
 */
template <typename Func> double gaussLegendreQuad(Func f, const double low, const double high)
{
    return boost::math::quadrature::gauss<double, 15>::integrate(f, low, high);
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
 *
 */
std::vector<std::vector<double>> choleskyDecomp(const std::vector<std::vector<double>> &matrix);

} // namespace util

#endif // UTIL_H
