#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <utility>
#include <vector>

#include "DecaySimulator.h"

#include "PullStudyHelpers.h"

/*
 * Test we find the correct values of a, b and c
 * Use a small value of tolerance; we are comparing resonably large floats
 */
BOOST_AUTO_TEST_CASE(test_expected_params, *boost::unit_test::tolerance(0.0000000001))
{
    // Set our decay parameters and what we expect a, b and c to evaluate to
    DecayParams_t DecayParams = DecayParams_t{.x = 1, .y = 2, .r = 3, .z_im = 4, .z_re = 5, .width = 6};

    // Parameters are described by eq. 2.10 in 1412.7254v2
    double expectedA = 9;
    double expectedB = 252;
    double expectedC = 45;

    std::vector<double> params = PullStudyHelpers::expectedParams(DecayParams);

    BOOST_CHECK(expectedA == params[0]);
    BOOST_CHECK(expectedB == params[1]);
    BOOST_CHECK(expectedC == params[2]);
}

/*
 * Test we find the correct relative number of DCS to CF decays
 */
BOOST_AUTO_TEST_CASE(test_numbers_of_decays)
{
    // TODO set some reasonable DecayParams_t that gives us a positive number of DCS/CF decays
    // Test that our code is right
}
