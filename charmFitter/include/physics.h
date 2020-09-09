/*
 * Functions representing physics stuff like decay rates, efficiencies
 */
#ifndef PHYSICS_H
#define PHYSICS_H

#include "fitterUtil.h"

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
 * Efficiency function
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
 * CF decay rate with efficiency, parametrised by tau
 */
inline double cfRate(const double time, const double width, const double tau)
{
    return exp(-1.0 * width * time) * efficiency(tau, time);
}

/*
 * CF decay rate with efficiency parametrised by tau
 */
inline double cfRate(const double time, const FitterUtil::DecayParams_t &decayParams, const double tau)
{
    return cfRate(time, decayParams.width, tau);
}

/*
 * Expected DCS decay rate at a given time with efficiency parameterise by tau
 */
inline double dcsRate(const double time, const FitterUtil::DecayParams_t &decayParams, const double tau)
{
    auto abcParams = expectedParams(decayParams);
    return rateRatio(time, abcParams) * cfRate(time, decayParams.width, tau);
}

/*
 * Expected DCS decay rate at a given time with efficiency parametrised by tau
 */
inline double dcsRate(const double time, const std::array<double, 3> &abcParams, const double width, const double tau)
{
    return rateRatio(time, abcParams) * cfRate(time, width, tau);
}

/*
 * Integral of CF rate between limits, with efficiency
 * Uses 1/width as tau in the efficiency
 *
 * errorEstimate can be set to get an estimate of the error in the integral
 */
inline double
cfIntegralWithEfficiency(const double low, const double high, const double width, const double efficiencyTimescale)
{
    // should really use std::bind
    auto f = [&](double x) { return cfRate(x, width, efficiencyTimescale); };
    return FitterUtil::gaussLegendreQuad(f, low, high);
}

/*
 * Integral of DCS rate between limits with efficiency
 * Uses 1/width as tau in the efficiency

 * errorEstimate can be set to get an estimate of the error in the integral
 */
inline double dcsIntegralWithEfficiency(const double                 low,
                                        const double                 high,
                                        const std::array<double, 3> &abcParams,
                                        const double                 width,
                                        const double                 efficiencyTimescale)
{
    // should really use std::bind
    auto f = [&](double x) { return dcsRate(x, abcParams, width, efficiencyTimescale); };
    return FitterUtil::gaussLegendreQuad(f, low, high);
}

/*
 * Find the number of DCS decays we need to simulate, given the number of CF decays and our phase space parameters.
 *
 * The ratio of DCS to CF decays is calculated from the ratio of the integrals of their decay rates.
 * CF rate is exponential; DCS is exp * (a + bt + ct^2)
 *
 * Returns a double so e.g. can be used as the mean of a distribution. Cast to an integer type before using as a
 * count!
 *
 */
double numDCSDecays(const size_t                     numCFDecays,
                    const FitterUtil::DecayParams_t &phaseSpaceParams,
                    double                           maxTime,
                    double                           efficiencyTimescale);

} // namespace Phys

#endif // PHYSICS_H
