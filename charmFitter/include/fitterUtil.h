#ifndef FITTER_UTIL_H
#define FITTER_UTIL_H

#include <vector>

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

} // namespace FitterUtil

#endif // FITTER_UTIL_H
