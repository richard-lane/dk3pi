#ifndef SIMULATOR_PYTHON_H
#define SIMULATOR_PYTHON_H

#include <utility>
#include <array>
#include <vector>

/*
 * Binding for expected params fcn
 */
std::array<double, 3> expectedParamsBinding(const std::array<double, 6>& decayParamsArr);

/*
 * Simulate decays, return them
 *
 * Generates n RS decays and the corresponding number of WS decays according to the provided decay params
 *
 * Decay params passed as an array {x, y, r, z_im, z_re, width}
 *
 */
std::pair<std::vector<double>, std::vector<double>>
simulate(const size_t n,
         const std::array<double, 6>& decayParamsArr,
         const double maxTime,
         const uint32_t seed);

#endif // SIMULATOR_PYTHON_H

