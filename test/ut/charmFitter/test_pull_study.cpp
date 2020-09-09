#include <boost/test/unit_test.hpp>

#include <cmath>
#include <utility>
#include <vector>

#include "DecaySimulator.h"
#include "fitterUtil.h"
#include "physics.h"

/*
 * Test we find the correct relative number of DCS to CF decays
 */
BOOST_AUTO_TEST_CASE(test_numbers_of_decays, *boost::unit_test::tolerance(1))
{
    // Set some reasonable DecayParams that gives us a positive number of DCS/CF decays
    // Gives a, b, c = 0.0025, 0.55, 101.5625
    FitterUtil::DecayParams_t DecayParams = {
        .x     = 0.004,
        .y     = 0.007,
        .r     = 0.05,
        .z_im  = -0.3,
        .z_re  = 0.8,
        .width = 2500.0,
    };

    // Test that our code is right
    size_t numCFDecays = 100000000;
    double maxTime     = 0.005;

    // Work this out using wolframalpha
    BOOST_CHECK(std::abs(Phys::numDCSDecays(numCFDecays, DecayParams, maxTime, 0) - 275248) < 1);
}
