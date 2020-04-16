#ifndef UTIL_H
#define UTIL_H

#include <boost/filesystem.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

#include "TObject.h"

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
 * Expected CF decay rate at a given time
 */
inline double rightSignDecayRate(const double time, const DecayParams_t &decayParams)
{
    return exp(-1.0 * decayParams.width * time);
}

/*
 * Expected DCS decay rate at a given time
 */
inline double wrongSignDecayRate(const double time, const DecayParams_t &decayParams)
{
    // Write the decay rate as (a + bt + ct^2)e^(-gamma*t) (ignoring overall factor of B^2 that has been taken out)
    double a = pow(decayParams.r, 2);
    double b =
        decayParams.r * (decayParams.y * decayParams.z_re + decayParams.x * decayParams.z_im) * decayParams.width;
    double c = 0.25 * (pow(decayParams.x, 2) + pow(decayParams.y, 2)) * pow(decayParams.width, 2);

    return (a + b * time + c * pow(time, 2)) * exp(-1.0 * decayParams.width * time);
}

/*
 * Integral of CF rate between limits
 */
inline double cfIntegral(const double low, const double high, const DecayParams_t &decayParams)
{
    // should really use std::bind
    auto f = [&](double x) { return rightSignDecayRate(x, decayParams); };
    return boost::math::quadrature::trapezoidal(f, low, high, 1e-10, 20);
}

/*
 * Integral of DCS rate between limits
 */
inline double dcsIntegral(const double low, const double high, const DecayParams_t &decayParams)
{
    // should really use std::bind
    auto f = [&](double x) { return wrongSignDecayRate(x, decayParams); };
    return boost::math::quadrature::trapezoidal(f, low, high, 1e-10, 20);
}

} // namespace util

#endif // UTIL_H
