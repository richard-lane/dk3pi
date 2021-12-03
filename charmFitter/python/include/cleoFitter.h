#ifndef CLEOFITTER_PYTHON_H
#define CLEOFITTER_PYTHON_H

#include <vector>
#include <array>

#include "CleoCombinationFitter.h"
#include "cleo_interface.h"

/*
 * Find CLEO likelihood at a series of Z values
 *
 * Return value indexed as (i + r * nImZVals) for imaginary+real (i, r)
 *
 */
std::vector<double> cleoLikelihoods(const std::vector<double>&   reZVals,
                                    const std::vector<double>&   imZVals,
                                    const std::array<double, 6>& decayParams,
                                    const int                    binNumber);

/*
 * Perform fits fixing Re(Z) and Im(Z) to the provided values; return an array of likelihoods
 *
 * Return value indexed as (i + r * nImZVals) for imaginary+real (i, r)
 *
 */
std::vector<double> cleoZScan(const std::vector<double>&  reZVals,
                             const std::vector<double>&   imZVals,
                             const std::vector<double>&   rsDecayTimes,
                             const std::vector<double>&   rsWeights,
                             const std::vector<double>&   wsDecayTimes,
                             const std::vector<double>&   wsWeights,
                             const std::vector<double>&   binLimits,
                             const std::array<double, 6>& initialVals,
                             const std::array<double, 6>& initialErrs,
                             const int                    binNumber);

#endif // CLEOFITTER_PYTHON_H
