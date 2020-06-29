/*
 * General utility functions and structs
 *
 * Should probably rename this
 */
#ifndef EFFICIENCY_UTIL_H
#define EFFICIENCY_UTIL_H

#include <boost/filesystem.hpp>

/*
 * Particle decay parameters
 */
typedef struct kinematicParams {
    double px{0};
    double py{0};
    double pz{0};
    double energy{0};
} kinematicParams_t;

/*
 * Parameters for decays D -> K pi1 pi2 pi3
 *
 * pi and pi2 should be the opposite charge of the K; p3 is the same charge as the K
 */
typedef struct dDecayParameters {
    kinematicParams_t dParams;

    kinematicParams_t kParams;
    kinematicParams_t pi1Params;
    kinematicParams_t pi2Params;
    kinematicParams_t pi3Params;

    /*
     * Tag the sign of the Kaon
     *
     * Defaults to K+
     */
    bool kPlus{true};
} dDecay_t;

#endif // EFFICIENCY_UTIL_H
