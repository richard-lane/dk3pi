#include "physics.h"
#include "util.h"

namespace Phys
{

double numDCSDecays(const size_t numCFDecays, const FitterUtil::DecayParams_t& phaseSpaceParams, double maxTime)
{
    // Lambdas for RS and WS rates
    auto rsRate = [&phaseSpaceParams](const double t) { return cfRate(t, phaseSpaceParams.width); };
    auto wsRate = [&phaseSpaceParams](const double t) { return dcsRate(t, phaseSpaceParams); };

    // Our formula is numDcs = numCf * (DCS integral / CF integral), where we integrate over all allowed times
    double dcsIntegral = util::gaussLegendreQuad(wsRate, 0, maxTime);
    double cfIntegral  = util::gaussLegendreQuad(rsRate, 0, maxTime);

    return numCFDecays * dcsIntegral / cfIntegral;
}

} // namespace Phys
