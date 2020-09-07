#include "fitterUtil.h"

namespace FitterUtil
{

std::vector<double> expectedParams(const DecayParams_t &phaseSpaceParams)
{
    double expected_a = phaseSpaceParams.r * phaseSpaceParams.r;
    double expected_b = phaseSpaceParams.r *
                        (phaseSpaceParams.y * phaseSpaceParams.z_re + phaseSpaceParams.x * phaseSpaceParams.z_im) *
                        phaseSpaceParams.width;
    double expected_c = 0.25 * (std::pow(phaseSpaceParams.x, 2) + std::pow(phaseSpaceParams.y, 2)) *
                        std::pow(phaseSpaceParams.width, 2);

    return std::vector<double>{expected_a, expected_b, expected_c};
}

std::vector<double> exponentialBinLimits(const double maxTime, const double decayConstant, const size_t numBins)
{
    std::vector<double> binLimits{};
    for (size_t i = 0; i <= numBins; ++i) {
        double x = (double)i / numBins;
        double z = 1 - std::exp(-1 * decayConstant * maxTime);
        binLimits.push_back((-1 / decayConstant) * std::log(1 - z * x));
    }

    return binLimits;
}


} // namespace FitterUtil
