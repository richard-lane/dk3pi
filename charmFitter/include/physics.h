/*
 * Functions representing physics stuff like decay rates, efficiencies
 */
#ifndef PHYSICS_H
#define PHYSICS_H

#include "fitterUtil.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <float.h>

#define INTEGRAL_TOLERANCE (std::numeric_limits<double>::epsilon())
#define INTEGRAL_REFINEMENTS (25)

namespace Phys
{
/*
 * Find our expected a, b and c in (a + bt + ct^2) from a set of phase space parameters.
 *
 * Returns an array {a, b, c}
 */
inline std::array<double, 3> expectedParams(const FitterUtil::DecayParams_t &phaseSpaceParams)
{
    double expected_a = phaseSpaceParams.r * phaseSpaceParams.r;
    double expected_b = phaseSpaceParams.r *
                        (phaseSpaceParams.y * phaseSpaceParams.z_re + phaseSpaceParams.x * phaseSpaceParams.z_im) *
                        phaseSpaceParams.width;
    double expected_c = 0.25 * (std::pow(phaseSpaceParams.x, 2) + std::pow(phaseSpaceParams.y, 2)) *
                        std::pow(phaseSpaceParams.width, 2);

    return {expected_a, expected_b, expected_c};
}
/*
 * Ratio of rates given parameters a, b and c
 */
inline double rateRatio(const double time, const std::array<double, 3> &abcParams)
{
    return abcParams[0] + time * abcParams[1] + time * time * abcParams[2];
}

/*
 * Ratio of rates given decay parameters
 */
inline double rateRatio(const double time, const FitterUtil::DecayParams_t &decayParams)
{
    return rateRatio(time, expectedParams(decayParams));
}

/*
 * Example efficiency function
 *
 * tanh(time/tau)
 */
inline double efficiency(const double tau, const double time)
{
    // Special case for tau = 0.0; efficiency is just 1 here
    if (std::abs(tau) <= DBL_EPSILON) {
        return 1.0;
    }
    return std::tanh(time / tau);
}

/*
 * CF decay rate
 */
inline double cfRate(const double time, const double width)
{
    return exp(-1.0 * width * time);
}

/*
 * CF decay rate with efficiency parametrised by tau
 */
inline double cfRate(const double time, const FitterUtil::DecayParams_t &decayParams)
{
    return cfRate(time, decayParams.width);
}

/*
 * Expected DCS decay rate at a given time with efficiency parameterise by tau
 */
inline double dcsRate(const double time, const FitterUtil::DecayParams_t &decayParams)
{
    auto abcParams = expectedParams(decayParams);
    return rateRatio(time, abcParams) * cfRate(time, decayParams.width);
}

/*
 * Expected DCS decay rate at a given time
 */
inline double dcsRate(const double time, const std::array<double, 3> &abcParams, const double width)
{
    return rateRatio(time, abcParams) * cfRate(time, width);
}

/*
 * DCS decay rate
 */
inline double dcsRate(const double time, const std::vector<double> &params)
{
    // Assuming parameter ordering here
    return dcsRate(time, FitterUtil::DecayParams_t{params[0], params[1], params[2], params[3], params[4], params[5]});
}

/*
 * DCS decay rate using polynomial parameters
 */
inline double dcsRatePoly(const double time, const std::vector<double> &params, const double width)
{
    std::array<double, 3> paramArray;
    std::copy_n(params.begin(), 3, paramArray.begin());
    return dcsRate(time, paramArray, width);
} // namespace Phys

/*
 * Find the number of DCS decays we need to simulate, given the number of CF decays and our phase space parameters.
 *
 * The ratio of DCS to CF decays is calculated from the ratio of the integrals of their decay rates.
 * CF rate is exponential; DCS is exp * (a + bt + ct^2)
 *
 * This function doesn't take efficiencies into account, since the ratio of RS/WS counts is agnostic to the form of the
 * LHCb efficiency
 *
 * Returns a double so e.g. can be used as the mean of a distribution. Cast to an integer type before using as a
 * count!
 *
 */
double numDCSDecays(const size_t numCFDecays, const FitterUtil::DecayParams_t &phaseSpaceParams, double maxTime);

} // namespace Phys

#endif // PHYSICS_H
