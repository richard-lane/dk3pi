#include "physics.h"

namespace Phys
{

double numDCSDecays(const size_t                     numCFDecays,
                    const FitterUtil::DecayParams_t &phaseSpaceParams,
                    double                           maxTime,
                    double                           efficiencyTimescale)
{
    // Our formula is numDcs = numCf * (DCS integral / CF integral), where we integrate over all allowed times
    double dcsIntegral = dcsIntegralWithEfficiency(
        0, maxTime, expectedParams(phaseSpaceParams), phaseSpaceParams.width, efficiencyTimescale);
    double cfIntegral = cfIntegralWithEfficiency(0, maxTime, phaseSpaceParams.width, efficiencyTimescale);

    return numCFDecays * dcsIntegral / cfIntegral;
}

} // namespace Phys
