#ifndef TEST_EFFICIENCY_UTIL_CPP
#define TEST_EFFICIENCY_UTIL_CPP

#define BOOST_TEST_DYN_LINK
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

#endif // TEST_EFFICIENCY_UTIL_CPP
