
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
 * Check that accept-reject works for a point on the RS curve
 * It's impossible to have a point rejected from this distribution
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

    SimulatedDecays MyDecays = SimulatedDecays(0.002, DecayParams);

    BOOST_CHECK(MyDecays.isAccepted(0.0005, 0.25, true));
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

    SimulatedDecays MyDecays = SimulatedDecays(0.002, DecayParams);

    BOOST_CHECK(MyDecays.isAccepted(0.0005, 0.699, false));
    BOOST_CHECK(!MyDecays.isAccepted(0.0005, 0.700, false));
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

    SimulatedDecays MyDecays = SimulatedDecays(0.002, DecayParams);

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

    SimulatedDecays MyDecays = SimulatedDecays(0.002, DecayParams);
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

    SimulatedDecays MyDecays = SimulatedDecays(0.002, DecayParams);
    MyDecays.findDcsDecayTimes(1);

    BOOST_CHECK_THROW(MyDecays.plotRates(std::vector<double>{0, 1}), D2K3PiException);
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

    SimulatedDecays MyDecays = SimulatedDecays(0.002, DecayParams);
    MyDecays.findCfDecayTimes(1);
    MyDecays.findDcsDecayTimes(1);

    BOOST_CHECK_THROW(MyDecays.plotRates(std::vector<double>{0.005, 0.002}), D2K3PiException);
}

#endif // TEST_SIMULATOR_CPP
