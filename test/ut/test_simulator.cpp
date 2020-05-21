
#ifndef TEST_SIMULATOR_CPP
#define TEST_SIMULATOR_CPP

#define BOOST_TEST_DYN_LINK
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <stdio.h>
#include <utility>
#include <vector>

#include "../pull_study/DecaySimulator.h"
#include "D2K3PiError.h"
#include "physics.h"

/*
 * Floating point tolerance for this UT module
 */
#define SIMULATOR_CPP_TOLERANCE 0.0000001

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

    SimulatedDecays MyDecays = SimulatedDecays(0.002, DecayParams, 0);

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

    SimulatedDecays MyDecays = SimulatedDecays(0.002, DecayParams, 0);
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

    SimulatedDecays MyDecays = SimulatedDecays(0.002, DecayParams, 0);
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

    SimulatedDecays MyDecays = SimulatedDecays(0.002, DecayParams, 0);
    MyDecays.findCfDecayTimes(1);
    MyDecays.findDcsDecayTimes(1);

    BOOST_CHECK_THROW(MyDecays.plotRates(std::vector<double>{0.005, 0.002}), D2K3PiException);
}

#endif // TEST_SIMULATOR_CPP
