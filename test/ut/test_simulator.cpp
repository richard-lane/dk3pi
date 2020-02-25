
#ifndef TEST_SIMULATOR_CPP
#define TEST_SIMULATOR_CPP

#define BOOST_TEST_DYN_LINK
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <stdio.h>
#include <utility>
#include <vector>

#include "D2K3PiError.h"
#include "DecaySimulator.h"

/*
 * Floating point tolerance for this UT module
 */
#define SIMULATOR_CPP_TOLERANCE 0.0000001

/*
 * Check that we get the right RS decay rate
 */
BOOST_AUTO_TEST_CASE(test_rs_rate, *boost::unit_test::tolerance(SIMULATOR_CPP_TOLERANCE))
{
    DecayParams_t DecayParams = {
        .x     = 0.004,
        .y     = 0.007,
        .r     = 0.05,
        .z_im  = -0.3,
        .z_re  = 0.8,
        .width = 2500.0,
    };

    // Define what we expect our RS decay rate to be
    auto expectedRSRate = [](const DecayParams_t &MyDecayParams, const double time) {
        return std::exp(-1 * MyDecayParams.width * time);
    };

    // Create a decay simulator, check our RS rate is right
    std::pair<double, double> allowedTimes = std::make_pair(0, 0.5);
    std::pair<double, double> allowedRates = std::make_pair(0, 1);

    SimulatedDecays MyDecays = SimulatedDecays(allowedTimes, allowedRates, DecayParams);

    BOOST_CHECK(expectedRSRate(DecayParams, 0.1) == MyDecays._rightSignDecayRate(0.1));
}

/*
 * Check that we get the right WS decay rate
 */
BOOST_AUTO_TEST_CASE(test_ws_rate, *boost::unit_test::tolerance(SIMULATOR_CPP_TOLERANCE))
{
    DecayParams_t DecayParams = {
        .x     = 0.004,
        .y     = 0.007,
        .r     = 0.05,
        .z_im  = -0.3,
        .z_re  = 0.8,
        .width = 2500.0,
    };

    // Define what we expect our RS decay rate to be
    auto expectedWSRate = [](const DecayParams_t &MyDecayParams, const double time) {
        return std::exp(-1 * MyDecayParams.width * time) *
               (MyDecayParams.r * MyDecayParams.r +
                MyDecayParams.width * MyDecayParams.r *
                    (MyDecayParams.y * MyDecayParams.z_re + MyDecayParams.x * MyDecayParams.z_im) * time +
                0.25 * (MyDecayParams.x * MyDecayParams.x + MyDecayParams.y * MyDecayParams.y) * MyDecayParams.width *
                    MyDecayParams.width * time * time);
    };

    // Create a decay simulator, check our RS rate is right
    std::pair<double, double> allowedTimes = std::make_pair(0, 0.5);
    std::pair<double, double> allowedRates = std::make_pair(0, 1);

    SimulatedDecays MyDecays = SimulatedDecays(allowedTimes, allowedRates, DecayParams);

    BOOST_CHECK(expectedWSRate(DecayParams, 0.01) == MyDecays._wrongSignDecayRate(0.01));
}

/*
 * Check that accept-reject works for two points on the RS curve
 */
BOOST_AUTO_TEST_CASE(test_acc_rej_rs)
{
    DecayParams_t DecayParams = {
        .x     = 0.004,
        .y     = 0.007,
        .r     = 0.05,
        .z_im  = -0.3,
        .z_re  = 0.8,
        .width = 2500.0,
    };

    std::pair<double, double> allowedTimes = std::make_pair(0, 0.002);
    std::pair<double, double> allowedRates = std::make_pair(0, 1.3);
    SimulatedDecays           MyDecays     = SimulatedDecays(allowedTimes, allowedRates, DecayParams);

    BOOST_CHECK(MyDecays.isAccepted(0.0005, 0.25, true));
    BOOST_CHECK(!MyDecays.isAccepted(0.0005, 0.35, true));
}

/*
 * Check that accept-reject works for two points on the RS curve
 * Uses the same two points as the RS test as the curves are basically the same (right..?)
 */
BOOST_AUTO_TEST_CASE(test_acc_rej_ws)
{
    DecayParams_t DecayParams = {
        .x     = 0.004,
        .y     = 0.007,
        .r     = 0.05,
        .z_im  = -0.3,
        .z_re  = 0.8,
        .width = 2500.0,
    };

    std::pair<double, double> allowedTimes = std::make_pair(0, 0.002);
    std::pair<double, double> allowedRates = std::make_pair(0, 1.3);
    SimulatedDecays           MyDecays     = SimulatedDecays(allowedTimes, allowedRates, DecayParams);

    BOOST_CHECK(MyDecays.isAccepted(0.0005, 0.0006, false));
    BOOST_CHECK(!MyDecays.isAccepted(0.0005, 0.0010, false));
}

/*
 * Test that attempting to plot the hists without setting WS or RS causes an error
 */
BOOST_AUTO_TEST_CASE(test_hist_error_ws_rs_not_set)
{
    DecayParams_t DecayParams = {
        .x     = 0.004,
        .y     = 0.007,
        .r     = 0.05,
        .z_im  = -0.3,
        .z_re  = 0.8,
        .width = 2500.0,
    };

    std::pair<double, double> allowedTimes = std::make_pair(0, 0.002);
    std::pair<double, double> allowedRates = std::make_pair(0, 1.3);
    SimulatedDecays           MyDecays     = SimulatedDecays(allowedTimes, allowedRates, DecayParams);

    BOOST_CHECK_THROW(MyDecays.plotRates(std::vector<double>{0, 1}), D2K3PiException);
}

/*
 * Test that attempting to plot the hists without setting RS causes an error
 */
BOOST_AUTO_TEST_CASE(test_hist_error_rs_not_set)
{
    DecayParams_t DecayParams = {
        .x     = 0.004,
        .y     = 0.007,
        .r     = 0.05,
        .z_im  = -0.3,
        .z_re  = 0.8,
        .width = 2500.0,
    };

    std::pair<double, double> allowedTimes = std::make_pair(0, 0.002);
    std::pair<double, double> allowedRates = std::make_pair(0, 1.3);
    SimulatedDecays           MyDecays     = SimulatedDecays(allowedTimes, allowedRates, DecayParams);
    MyDecays.findCfDecayTimes(1);

    BOOST_CHECK_THROW(MyDecays.plotRates(std::vector<double>{0, 1}), D2K3PiException);
}

/*
 * Test that attempting to plot the hists without setting WS causes an error
 */
BOOST_AUTO_TEST_CASE(test_hist_error_ws_not_set)
{
    DecayParams_t DecayParams = {
        .x     = 0.004,
        .y     = 0.007,
        .r     = 0.05,
        .z_im  = -0.3,
        .z_re  = 0.8,
        .width = 2500.0,
    };

    std::pair<double, double> allowedTimes = std::make_pair(0, 0.002);
    std::pair<double, double> allowedRates = std::make_pair(0, 1.3);
    SimulatedDecays           MyDecays     = SimulatedDecays(allowedTimes, allowedRates, DecayParams);
    MyDecays.findDcsDecayTimes(1);

    BOOST_CHECK_THROW(MyDecays.plotRates(std::vector<double>{0, 1}), D2K3PiException);
}

/*
 * Test that using bin limits that dont cover the full range of times cause an error
 */
BOOST_AUTO_TEST_CASE(test_hist_bad_bin_limits)
{
    DecayParams_t DecayParams = {
        .x     = 0.004,
        .y     = 0.007,
        .r     = 0.05,
        .z_im  = -0.3,
        .z_re  = 0.8,
        .width = 2500.0,
    };

    std::pair<double, double> allowedTimes = std::make_pair(0, 0.002);
    std::pair<double, double> allowedRates = std::make_pair(0, 1.3);
    SimulatedDecays           MyDecays     = SimulatedDecays(allowedTimes, allowedRates, DecayParams);
    MyDecays.findCfDecayTimes(1);
    MyDecays.findDcsDecayTimes(1);

    BOOST_CHECK_THROW(MyDecays.plotRates(std::vector<double>{0, 0.0015}), D2K3PiException);
    BOOST_CHECK_THROW(MyDecays.plotRates(std::vector<double>{0.0005, 0.002}), D2K3PiException);
}

/*
 * Test that using unsorted bin limits causes an error
 */
BOOST_AUTO_TEST_CASE(test_hist_unsorted_bin_limits)
{
    DecayParams_t DecayParams = {
        .x     = 0.004,
        .y     = 0.007,
        .r     = 0.05,
        .z_im  = -0.3,
        .z_re  = 0.8,
        .width = 2500.0,
    };

    std::pair<double, double> allowedTimes = std::make_pair(0, 0.002);
    std::pair<double, double> allowedRates = std::make_pair(0, 1.3);
    SimulatedDecays           MyDecays     = SimulatedDecays(allowedTimes, allowedRates, DecayParams);
    MyDecays.findCfDecayTimes(1);
    MyDecays.findDcsDecayTimes(1);

    BOOST_CHECK_THROW(MyDecays.plotRates(std::vector<double>{0.005, 0.002}), D2K3PiException);
}

#endif // TEST_SIMULATOR_CPP
