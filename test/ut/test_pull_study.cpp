#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <utility>
#include <vector>

#include "../pull_study/DecaySimulator.h"
#include "../pull_study/PullStudyHelpers.h"

/*
 * Test mean + std dev function
 */
BOOST_AUTO_TEST_CASE(test_mean_std)
{
    std::vector<double>       data{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    std::pair<double, double> expected = std::make_pair(6.5, 3.45205);
    std::pair<double, double> actual   = PullStudyHelpers::meanAndStdDev(data);

    BOOST_CHECK_SMALL(expected.first - actual.first, 1e-5);
    BOOST_CHECK_SMALL(expected.second - actual.second, 1e-5);
}

/*
 * Test we find the correct relative number of DCS to CF decays
 */
BOOST_AUTO_TEST_CASE(test_numbers_of_decays, *boost::unit_test::tolerance(1))
{
    // Set some reasonable DecayParams that gives us a positive number of DCS/CF decays
    // Gives a, b, c = 0.0025, 0.55, 101.5625
    DecayParams_t DecayParams = {
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
    BOOST_CHECK(std::abs(PullStudyHelpers::numDCSDecays(numCFDecays, DecayParams, maxTime, 0) - 275248) < 1);
}
