#include "efficiencyUtil.h"

#include <vector>

#include <TGenPhaseSpace.h>

/*
 * Generate an event final state randomly, uniformly distributed over phase space
 *
 * The TGenPhaseSpace should have already been told what particle is decaying/what the decay products are; this function
 * doesn't check that this has been done!
 *
 * Need to provide a random number generator that generates numbers between 0 and the maximum weight of the phase space.
 *
 * Optionally returns the weight of the event generated as an OUT argument cus I'm really a C programmer
 */
std::vector<kinematicParams_t> randomEvent(TGenPhaseSpace&                         phaseSpace,
                                           std::mt19937* const                     generator,
                                           std::uniform_real_distribution<double>& uniformDistribution,
                                           double*                                 weight = nullptr);

/*
 * Generate a vector of four-body D->K3pi events, uniformly distributed in phase space
 *
 * Boolean flag indicates whether the kaon is positive
 */
std::vector<dDecay_t> flatDk3pi(const size_t numEvents, std::mt19937* const generator, const bool kPlus = true);
