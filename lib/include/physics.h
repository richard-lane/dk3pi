/*
 * Functions representing physics stuff like decay rates, efficiencies
 */
#ifndef PHYSICS_H
#define PHYSICS_H

#include <vector>

#include "util.h"

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
    std::vector<double> abcParams = util::expectedParams(decayParams);
    return rateRatio(time, abcParams) * rightSignDecayRate(time, decayParams.width);
}

/*
 * Expected DCS decay rate at a given time
 */
inline double wrongSignDecayRate(const double time, const std::vector<double> &abcParams, const double width)
{
    return rateRatio(time, abcParams) * rightSignDecayRate(time, width);
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
inline double cfRateWithEfficiency(const double time, const double width, const double tau)
{
    return rightSignDecayRate(time, width) * efficiency(tau, time);
}

/*
 * CF decay rate with efficiency parametrised by tau
 */
inline double cfRateWithEfficiency(const double time, const DecayParams_t &decayParams, const double tau)
{
    return cfRateWithEfficiency(time, decayParams.width, tau);
}

/*
 * Expected DCS decay rate at a given time with efficiency parameterise by tau
 */
inline double dcsRateWithEfficiency(const double time, const DecayParams_t &decayParams, const double tau)
{
    return wrongSignDecayRate(time, decayParams) * efficiency(tau, time);
}

/*
 * Expected DCS decay rate at a given time with efficiency parametrised by tau
 */
inline double
dcsRateWithEfficiency(const double time, const std::vector<double> &abcParams, const double width, const double tau)
{
    return wrongSignDecayRate(time, abcParams, width) * efficiency(tau, time);
}

/*
 * Integral of CF rate between limits, with efficiency
 * Uses 1/width as tau in the efficiency
 *
 * errorEstimate can be set to get an estimate of the error in the integral
 */
inline double cfIntegralWithEfficiency(const double low,
                                       const double high,
                                       const double width,
                                       const double efficiencyTimescale,
                                       const double tolerance      = INTEGRAL_TOLERANCE,
                                       const size_t maxRefinements = INTEGRAL_REFINEMENTS,
                                       double *     errorEstimate  = nullptr)
{
    // should really use std::bind
    auto f = [&](double x) { return cfRateWithEfficiency(x, width, efficiencyTimescale); };
    return boost::math::quadrature::trapezoidal(f, low, high, tolerance, maxRefinements, errorEstimate);
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
                                        const double               efficiencyTimescale,
                                        const double               tolerance      = INTEGRAL_TOLERANCE,
                                        const size_t               maxRefinements = INTEGRAL_REFINEMENTS,
                                        double *                   errorEstimate  = nullptr)
{
    // should really use std::bind
    auto f = [&](double x) { return dcsRateWithEfficiency(x, abcParams, width, efficiencyTimescale); };
    return boost::math::quadrature::trapezoidal(f, low, high, tolerance, maxRefinements, errorEstimate);
}

} // namespace Phys

#endif // PHYSICS_H
