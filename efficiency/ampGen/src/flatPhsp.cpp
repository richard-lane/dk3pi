#include "flatPhsp.h"
#include "efficiencyUtil.h"

#include <cassert>

#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>

std::vector<kinematicParams_t> randomEvent(TGenPhaseSpace&                         phaseSpace,
                                           std::mt19937* const                     generator,
                                           std::uniform_real_distribution<double>& uniformDistribution,
                                           double*                                 weight)
{
    double maxWeight   = phaseSpace.GetWtMax();
    double rndNumber   = 0;
    double eventWeight = 0; // These get re-set immediately

    // Vector for storing the kinematics of our decay particles
    size_t                         numFinalStateParticles = phaseSpace.GetNt();
    std::vector<kinematicParams_t> eventKinematics(numFinalStateParticles);

    do {
        // Generate a random number (between 0 and maxWeight) and an event
        rndNumber   = uniformDistribution(*generator);
        eventWeight = phaseSpace.Generate();
    } while (rndNumber > eventWeight);

    // Try and catch a potential known ROOT bug
    assert(rndNumber < maxWeight);
    if (eventWeight > maxWeight) {
        throw TGenPhspBug("Weight generated > max weight expected");
    }

    // Fill in our vector of decay kinematics
    for (size_t i = 0; i < numFinalStateParticles; ++i) {
        TLorentzVector particleData = *phaseSpace.GetDecay(i);
        eventKinematics[i].px       = particleData.Px();
        eventKinematics[i].py       = particleData.Py();
        eventKinematics[i].pz       = particleData.Pz();
        eventKinematics[i].energy   = particleData.E();
    }

    // Set weight
    if (weight) {
        *weight = eventWeight;
    }

    return eventKinematics;
}

std::vector<dDecay_t> flatDk3pi(const size_t numEvents, std::mt19937* const generator, const bool kPlus)
{
    // Initialise vector to the right length
    std::vector<dDecay_t> flatEvents(numEvents);

    // Our decaying particle should be a D meson at rest and the decay products are K 3pi
    TLorentzVector        dMomentum(0.0, 0.0, 0.0, D_MASS_GEV);
    std::array<double, 4> finalStateMasses = {K_MASS_GEV, PI_MASS_GEV, PI_MASS_GEV, PI_MASS_GEV};

    // Vector for storing our kinematic parameters (as this is what randomEvent returns)
    std::vector<kinematicParams_t> kinematicParams(4);

    TGenPhaseSpace PhaseSpace;
    PhaseSpace.SetDecay(dMomentum, finalStateMasses.size(), finalStateMasses.data());

    // We will need to generate random numbers between 0 and the max weight of our phase space
    double                                 maxWeight = PhaseSpace.GetWtMax();
    std::uniform_real_distribution<double> uniformDistribution(0.0, maxWeight);

    for (size_t i = 0; i < numEvents; ++i) {
        kinematicParams = randomEvent(PhaseSpace, generator, uniformDistribution);

        // The other D params should be correctly default initialised to 0
        flatEvents[i].dParams.energy = dMomentum.E();

        flatEvents[i].kParams   = kinematicParams[0];
        flatEvents[i].pi1Params = kinematicParams[1];
        flatEvents[i].pi2Params = kinematicParams[2];
        flatEvents[i].pi3Params = kinematicParams[3];

        flatEvents[i].kPlus = kPlus;
    }

    return flatEvents;
}
