
/*
 * Mostly copy and pasted from the BES likelihood
 */

#ifndef CHARM_INTERFACE_H
#define CHARM_INTERFACE_H

#include "bes_interface.h"
#include "cleo_interface.h"
#include "fitterUtil.h"

#include <map>

namespace charm_threshhold
{

/*
 * Find the combined CLEO + BES-III chi2

 */
inline double chi2(const BES::Bin phspBin, const FitterUtil::DecayParams_t& params)
{
    const auto besChi2{BES::besLikelihood(phspBin, params)};

    const auto cleoLikelihood{CLEO::cleoLikelihood(static_cast<CLEO::Bin>(phspBin), params)};
    const auto cleoChi2{-2.0 * cleoLikelihood};

    // std::cout << cleoChi2 << " " << besChi2 << std::endl;

    return cleoChi2 + besChi2;
}

} // namespace charm_threshhold

#endif // CHARM_INTERFACE_H
