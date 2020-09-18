#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <random>

#include <TFile.h>
#include <TGenPhaseSpace.h>
#include <TH2D.h>

#include "ReadRoot.h"
#include "amplitudes.h"
#include "efficiencyUtil.h"
#include "flatPhsp.h"
#include "util.h"

BOOST_AUTO_TEST_SUITE(test_phase_space)

/*
 * Check that a randomly generated 2-body decay (flat in phsp) is produced with the correct invariant mass
 */
BOOST_AUTO_TEST_CASE(test_two_body_decay)
{
    // Decay is D -> Kpi
    TLorentzVector        dMomentum(0.0, 0.0, 0.0, D_MASS_GEV);
    std::array<double, 2> finalStateMasses = {K_MASS_GEV, PI_MASS_GEV};
    TGenPhaseSpace        PhaseSpace;
    PhaseSpace.SetDecay(dMomentum, finalStateMasses.size(), finalStateMasses.data());

    // Random number stuff
    std::random_device                     rd;
    std::mt19937                           gen(rd());
    std::uniform_real_distribution<double> uniformDistribution(0, PhaseSpace.GetWtMax());

    // Generate 5 events, check all of them have the right invariant mass
    for (size_t i = 0; i < 5; ++i) {
        std::vector<kinematicParams_t> event = randomEvent(PhaseSpace, &gen, uniformDistribution);
        BOOST_CHECK(invariantMass(event) == D_MASS_GEV);
    }
}

BOOST_AUTO_TEST_SUITE_END()
