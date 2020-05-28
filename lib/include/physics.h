/*
 * Functions representing physics stuff like decay rates, efficiencies
 */
#ifndef PHYSICS_H
#define PHYSICS_H

#include <vector>

#include "util.h"

#define INTEGRAL_TOLERANCE (std::numeric_limits<double>::epsilon())
#define INTEGRAL_REFINEMENTS (25)

namespace Phys
{
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
    return rateRatio(time, util::expectedParams(decayParams));
}

/*
 * Efficiency function
 *
 * tanh(time/tau)
 */
inline double efficiency(const double tau, const double time)
{
    // Special case for tau = 0.0; efficiency is just 1 here
    if (tau == 0.0) {
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
inline double cfRate(const double time, const DecayParams_t &decayParams, const double tau)
{
    return cfRate(time, decayParams.width, tau);
}

/*
 * Expected DCS decay rate at a given time with efficiency parameterise by tau
 */
inline double dcsRate(const double time, const DecayParams_t &decayParams, const double tau)
{
    std::vector<double> abcParams = util::expectedParams(decayParams);
    return rateRatio(time, abcParams) * cfRate(time, decayParams.width, tau);
}

/*
 * Expected DCS decay rate at a given time with efficiency parametrised by tau
 */
inline double dcsRate(const double time, const std::vector<double> &abcParams, const double width, const double tau)
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
    return util::gaussLegendreQuad(f, low, high);
}

/*
 * Integral of DCS rate between limits with efficiency
 * Uses 1/width as tau in the efficiency

 * errorEstimate can be set to get an estimate of the error in the integral
 */
inline double dcsIntegralWithEfficiency(const double               low,
                                        const double               high,
                                        const std::vector<double> &abcParams,
                                        const double               width,
                                        const double               efficiencyTimescale)
{
    // should really use std::bind
    auto f = [&](double x) { return dcsRate(x, abcParams, width, efficiencyTimescale); };
    return util::gaussLegendreQuad(f, low, high);
}

} // namespace Phys

#endif // PHYSICS_H
