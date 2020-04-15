#include <random>
#include <vector>

#include "testutil.h"

std::vector<double>
idealRatios(const std::vector<double>& times, const double error, const double a, const double b, const double c)
{
    std::random_device               randomEngine;
    std::mt19937                     mt(randomEngine());
    std::normal_distribution<double> distribution(0.0, 1.0);

    std::vector<double> ratios      = std::vector<double>(times.size(), -1);
    std::vector<double> ratioErrors = std::vector<double>(times.size(), -1);

    // Create idealised plot
    for (size_t i = 0; i < times.size(); ++i) {
        ratios[i] = ratio(a, b, c, times[i]) * (1 + error * distribution(mt));
    }
    return ratios;
}
