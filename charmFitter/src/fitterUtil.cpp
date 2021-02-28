#include "fitterUtil.h"

#include <cmath>

namespace FitterUtil
{

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
