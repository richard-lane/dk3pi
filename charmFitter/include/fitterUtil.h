#ifndef FITTER_UTIL_H
#define FITTER_UTIL_H

#include <vector>

#include <boost/math/quadrature/gauss.hpp>

#include <TObject.h>

namespace FitterUtil
{

/*
 * Struct encapsulating the parameters needed to simulate decays.
 * Mixing params x, y
 * Phase space params r, Im(Z) and Re(Z)
 * Particle data \Gamma for the D meson
 *
 * Parameters are all initalised to zero by default and should be set by the user.
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
 * Use Gauss-Legendre quadrature to find an approximation to the integral of f between low and high limits
 *
 * Evaluates the function at 15 points, as the weights and abcissa have been precalculated for this number of points.
 */
template <typename Func> double gaussLegendreQuad(Func f, const double low, const double high)
{
    return boost::math::quadrature::gauss<double, 15>::integrate(f, low, high);
}

} // namespace FitterUtil

#endif // FITTER_UTIL_H
