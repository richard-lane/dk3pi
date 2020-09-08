#include <boost/test/unit_test.hpp>

#include "flatPhsp.h"

/*
 * Test that if specified, our weight get set
 */
BOOST_AUTO_TEST_CASE(test_weight)
{
    // 2 body phase space
    std::array<double, 2> finalStateMasses = {0.3, 0.3};
    TGenPhaseSpace        PhaseSpace;
    TLorentzVector        initialMomm(0.0, 0.0, 0.0, 1.0);
    PhaseSpace.SetDecay(initialMomm, finalStateMasses.size(), finalStateMasses.data());

    // Random number generator + distribution
    std::random_device                     rd;
    std::mt19937                           gen(rd());
    std::uniform_real_distribution<double> uniformDistribution(0.0, PhaseSpace.GetWtMax());

    // Weight that could be set by the random event generation function
    double weight = 0;

    BOOST_CHECK_NO_THROW(randomEvent(PhaseSpace, &gen, uniformDistribution));
    BOOST_CHECK_NO_THROW(randomEvent(PhaseSpace, &gen, uniformDistribution, &weight));

    // Check that our weight has been set
    BOOST_CHECK(weight > 0);
}
