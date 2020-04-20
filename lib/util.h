#ifndef UTIL_H
#define UTIL_H

#include <boost/filesystem.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

#include "D2K3PiError.h"

#include "TObject.h"

#define INTEGRAL_TOLERANCE 1e-15L
#define INTEGRAL_REFINEMENTS 20

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
 * Combine a directory, filename and extension into a single boost:filesystem::path object.
 *
 * Effectively just does plotDir + plotName + fileExtension in an OS-intelligent way.
 *
 */
boost::filesystem::path concatPaths(std::string plotDir, std::string plotName, std::string fileExtension);

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
 * TODO: this currently sucks for the above reason. so make it not suck
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
 * Find our expected a, b and c in (a + bt + ct^2) from a set of phase space parameters.
 *
 * Returns a vector of {a, b, c}
 */
std::vector<double> expectedParams(const DecayParams_t &phaseSpaceParams);

/*
 * Ratio of rates given parameters a, b and c
 */
inline double rateRatio(const double time, const std::vector<double> &abcParams)
{
    return abcParams[0] + time * abcParams[1] + time * time * abcParams[2];
}

/*
 * Ratio of rates given decay parameters
 */
inline double rateRatio(const double time, const DecayParams_t &decayParams)
{
    return rateRatio(time, expectedParams(decayParams));
}

/*
 * Expected CF decay rate at a given time
 */
inline double rightSignDecayRate(const double time, const double width)
{
    return exp(-1.0 * width * time);
}

/*
 * Expected CF decay rate at a given time
 */
inline double rightSignDecayRate(const double time, const DecayParams_t &decayParams)
{
    return rightSignDecayRate(time, decayParams.width);
}

/*
 * Expected DCS decay rate at a given time
 */
inline double wrongSignDecayRate(const double time, const DecayParams_t &decayParams)
{
    // Write the decay rate as (a + bt + ct^2)e^(-gamma*t) (ignoring overall factor of B^2 that has been taken out)
    std::vector<double> abcParams = expectedParams(decayParams);

    return (abcParams[0] + abcParams[1] * time + abcParams[2] * time * time) * exp(-1.0 * decayParams.width * time);
}

/*
 * Integral of CF rate between limits
 *
 * errorEstimate can be set to get an estimate of the error in the integral
 */
inline double cfIntegral(const double low,
                         const double high,
                         const double width,
                         const double tolerance      = INTEGRAL_TOLERANCE,
                         const size_t maxRefinements = INTEGRAL_REFINEMENTS,
                         double *     errorEstimate  = nullptr)
{
    // should really use std::bind
    auto f = [&](double x) { return rightSignDecayRate(x, width); };

    return boost::math::quadrature::trapezoidal(f, low, high, tolerance, maxRefinements, errorEstimate);
}

/*
 * Integral of DCS rate between limits
 *
 * errorEstimate can be set to get an estimate of the error in the integral
 */
inline double dcsIntegral(const double               low,
                          const double               high,
                          const std::vector<double> &abcParams,
                          const double               width,
                          const double               tolerance      = INTEGRAL_TOLERANCE,
                          const size_t               maxRefinements = INTEGRAL_REFINEMENTS,
                          double *                   errorEstimate  = nullptr)
{
    // should really use std::bind
    auto f = [&](double x) { return (abcParams[0] + abcParams[1] * x + abcParams[2] * x * x) * exp(-1.0 * width * x); };

    return boost::math::quadrature::trapezoidal(f, low, high, tolerance, maxRefinements, errorEstimate);
}

/*
 * Integral of CF rate between limits
 *
 * errorEstimate can be set to get an estimate of the error in the integral
 */
inline double cfIntegral(const double         low,
                         const double         high,
                         const DecayParams_t &decayParams,
                         const double         tolerance      = INTEGRAL_TOLERANCE,
                         const size_t         maxRefinements = INTEGRAL_REFINEMENTS,
                         double *             errorEstimate  = nullptr)
{
    return cfIntegral(low, high, decayParams.width, tolerance, maxRefinements, errorEstimate);
}

/*
 * Integral of DCS rate between limits
 *
 * errorEstimate can be set to get an estimate of the error in the integral
 */
inline double dcsIntegral(const double         low,
                          const double         high,
                          const DecayParams_t &decayParams,
                          const double         tolerance      = INTEGRAL_TOLERANCE,
                          const size_t         maxRefinements = INTEGRAL_REFINEMENTS,
                          double *             errorEstimate  = nullptr)
{
    return dcsIntegral(
        low, high, util::expectedParams(decayParams), decayParams.width, tolerance, maxRefinements, errorEstimate);
}

/*
 * Analytical integral of CF rate between limits
 */
inline double analyticalCfIntegral(const double low, const double high, const DecayParams_t &decayParams)
{
    return (std::exp(-decayParams.width * low) - std::exp(-decayParams.width * high)) / decayParams.width;
}

/*
 * Analytical integral of DCS rate between limits
 */
inline double analyticalDcsIntegral(const double low, const double high, const DecayParams_t &decayParams)
{
    // Our integral is composed of three terms, X + Y + Z
    std::vector<double> abcParams = expectedParams(decayParams);

    double X = abcParams[0] * analyticalCfIntegral(low, high, decayParams);

    auto y = [&](double x) {
        return std::exp(-decayParams.width * x) * (decayParams.width * x + 1) / (decayParams.width * decayParams.width);
    };
    double Y = abcParams[1] * (y(low) - y(high));

    auto z = [&](double x) {
        return std::exp(-decayParams.width * x) *
               (2 + 2 * decayParams.width * x + decayParams.width * decayParams.width * x * x) /
               (decayParams.width * decayParams.width * decayParams.width);
    };
    double Z = abcParams[2] * (z(low) - z(high));

    return X + Y + Z;
}

} // namespace util

#endif // UTIL_H
