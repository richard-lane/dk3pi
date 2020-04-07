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

    std::vector<double> params = util::expectedParams(DecayParams);

    BOOST_CHECK(expectedA == params[0]);
    BOOST_CHECK(expectedB == params[1]);
    BOOST_CHECK(expectedC == params[2]);
}

/*
 * Test we find the correct relative number of DCS to CF decays
 */
BOOST_AUTO_TEST_CASE(test_numbers_of_decays)
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
    size_t              numCFDecays = 100000000;
    double              maxTime     = 0.005;
    std::vector<double> params      = util::expectedParams(DecayParams);
    double              a           = params[0];
    double              b           = params[1];
    double              c           = params[2];

    // alpha, beta, gamma are the three parts of the integral of (a + bt + ct2)exp(...)
    double alpha = (a / DecayParams.width) * (1 - std::exp(-1 * DecayParams.width * maxTime));
    double beta  = (b / (DecayParams.width * DecayParams.width)) *
                  (1 - (DecayParams.width * maxTime + 1) * std::exp(-1 * DecayParams.width * maxTime));
    double gamma = (c / std::pow(DecayParams.width, 3)) *
                   (2 - (2 + DecayParams.width * maxTime * (DecayParams.width * maxTime + 2))) *
                   std::exp(-1 * DecayParams.width * maxTime);

    double expectedNumDCS =
        numCFDecays * (alpha + beta + gamma) * DecayParams.width / (1 - std::exp(-1 * DecayParams.width * maxTime));

    // There might be rounding errors so check they're within 2 of each other... that should be good enough...
    BOOST_CHECK(std::abs(expectedNumDCS - PullStudyHelpers::numDCSDecays(numCFDecays, DecayParams, maxTime)) < 2.5);
}
