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
    // Write the decay rate as (a + bt + ct^2)e^(-gamma*t) (ignoring overall factor of B^2 that has been taken out)
    std::vector<double> abcParams = util::expectedParams(decayParams);

    return (abcParams[0] + abcParams[1] * time + abcParams[2] * time * time) * exp(-1.0 * decayParams.width * time);
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
 * Analytical integral of CF rate between limits
 */
inline double analyticalCfIntegral(const double low, const double high, const double width)
{
    return (std::exp(-width * low) - std::exp(-width * high)) / width;
}

/*
 * Analytical integral of DCS rate between limits
 */
inline double
analyticalDcsIntegral(const double low, const double high, const std::vector<double> &abcParams, const double width)
{
    // Our integral is composed of three terms, X + Y + Z

    double X = abcParams[0] * analyticalCfIntegral(low, high, width);

    auto   y = [&](double x) { return std::exp(-width * x) * (width * x + 1) / (width * width); };
    double Y = abcParams[1] * (y(low) - y(high));

    auto z = [&](double x) {
        return std::exp(-width * x) * (2 + 2 * width * x + width * width * x * x) / (width * width * width);
    };
    double Z = abcParams[2] * (z(low) - z(high));

    return X + Y + Z;
}

/*
 * Analytical integral of CF rate between limits
 */
inline double analyticalCfIntegral(const double low, const double high, const DecayParams_t &decayParams)
{
    return analyticalCfIntegral(low, high, decayParams.width);
}

/*
 * Analytical integral of DCS rate between limits
 */
inline double analyticalDcsIntegral(const double low, const double high, const DecayParams_t &decayParams)
{
    std::vector<double> abcParams = util::expectedParams(decayParams);
    return analyticalDcsIntegral(low, high, abcParams, decayParams.width);
}

/*
 * Efficiency function
 * Step function at tau
 */
inline double efficiency(const double tau, const double time)
{
    return time < tau ? 0.0 : 1.0;
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

} // namespace Phys

#endif // PHYSICS_H
