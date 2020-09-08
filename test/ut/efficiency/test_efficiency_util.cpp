#include <boost/test/unit_test.hpp>

#include "efficiencyUtil.h"

/*
 * Test invariant mass calculation for three bodies
 */
BOOST_AUTO_TEST_CASE(test_invariant_mass)
{
    kinematicParams_t particle1{.px = 1, .py = 2, .pz = 3, .energy = 10};
    kinematicParams_t particle2{.px = 0, .py = -1, .pz = -2, .energy = 5};
    kinematicParams_t particle3{.px = 1, .py = 1, .pz = 1, .energy = 7};

    double expectedInvariantMass = 21.72556098;

    BOOST_CHECK_CLOSE(
        invariantMass(std::vector<kinematicParams_t>{particle1, particle2, particle3}), expectedInvariantMass, 1e-7);
}

/*
 * Test calculation taking an event to a phsp location
 */
BOOST_AUTO_TEST_CASE(test_event_to_phsp)
{
    kinematicParams_t d{.px = 1, .py = 2, .pz = 3, .energy = 10};
    kinematicParams_t k{.px = 0, .py = -1, .pz = -2, .energy = 5};
    kinematicParams_t pi1{.px = 1, .py = 1, .pz = 1, .energy = 2};
    kinematicParams_t pi2{.px = 4, .py = 5, .pz = -2, .energy = 23};
    kinematicParams_t pi3{.px = 2, .py = 7, .pz = -6, .energy = 24};

    dDecay_t event{.dParams = d, .kParams = k, .pi1Params = pi1, .pi2Params = pi2, .pi3Params = pi3, .kPlus = true};

    std::vector<double> expectedInvariantMasses{invariantMass(std::vector<kinematicParams_t>{k, pi1}),
                                                invariantMass(std::vector<kinematicParams_t>{pi1, pi2}),
                                                invariantMass(std::vector<kinematicParams_t>{pi2, pi3}),
                                                invariantMass(std::vector<kinematicParams_t>{k, pi1, pi2}),
                                                invariantMass(std::vector<kinematicParams_t>{pi1, pi2, pi3})};
    std::vector<double> invariantMasses = event2invariantMasses(event);

    for (size_t i = 0; i < invariantMasses.size(); ++i) {
        BOOST_CHECK_CLOSE(invariantMasses[i], expectedInvariantMasses[i], 1e-7);
    }
}

/*
 * Test pT
 */
BOOST_AUTO_TEST_CASE(test_pt)
{
    kinematicParams_t particle1{.px = 0, .py = 1, .pz = 3, .energy = 10};
    kinematicParams_t particle2{.px = 2, .py = 3, .pz = 3, .energy = 10};
    double            expectedPT = std::sqrt(20);

    std::vector<kinematicParams_t> particles{particle1, particle2};

    BOOST_CHECK_CLOSE(expectedPT, pT(particles), 1e-7);
}
