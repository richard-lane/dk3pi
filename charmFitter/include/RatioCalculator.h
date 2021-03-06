/*
 * Bin decay times into specified time bins and calculate the ratio of CF to DCS decays in each.
 */
#ifndef RATIOCALCULATOR_H
#define RATIOCALCULATOR_H

#include <utility>
#include <vector>

#include "D2K3PiError.h"

namespace RatioCalculator
{
/*
 * Dataset contains zeros
 */
struct ZerosInData : public D2K3PiException {
    const char *what() const throw() { return "Zeros in dataset; cannot take ratio"; }
};

/*
 * Datasets of different sizes
 */
struct WrongSizeDataset : public D2K3PiException {
    const char *what() const throw() { return "Numerator and denominator have different numbers of points"; }
};

/*
 * Calculate the ratio of two binned datasets and the associated error, assuming they obey Poisson statistics
 *
 * Returns a pair <ratios, errors>
 */
std::pair<std::vector<double>, std::vector<double>> ratioAndError(const std::vector<size_t> &denominator,
                                                                  const std::vector<size_t> &numerator);
} // namespace RatioCalculator

#endif // RATIOCALCULATOR_H
