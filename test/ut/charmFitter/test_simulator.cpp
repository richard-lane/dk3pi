#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <memory>
#include <stdio.h>
#include <utility>
#include <vector>

#include "D2K3PiError.h"
#include "DecaySimulator.h"
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
    auto                          fcn = [](double x) { return 2 * x; };
    auto                          gen = [&]() { return 1.; };
    std::shared_ptr<std::mt19937> generator;

    SimulatedDecays MyDecays = SimulatedDecays(gen, fcn, fcn, fcn, std::make_pair(0.0, 1.), generator);

    BOOST_CHECK_THROW(MyDecays.plotRates(std::vector<double>{0, 1}), D2K3PiException);
}

/*
 * Test that using unsorted bin limits causes an error
 */
BOOST_AUTO_TEST_CASE(test_hist_unsorted_bin_limits)
{
    auto                          rate      = [](double x) { return 2 * x; };
    auto                          gen       = [&]() { return 1.; };
    std::shared_ptr<std::mt19937> generator = std::make_shared<std::mt19937>();

    SimulatedDecays MyDecays = SimulatedDecays(gen, rate, rate, rate, std::make_pair(0., 1.), generator);
    MyDecays.findCfDecayTimes(1);
    MyDecays.findDcsDecayTimes(1);

    BOOST_CHECK_THROW(MyDecays.plotRates(std::vector<double>{0.005, 0.002}), D2K3PiException);
}

/*
 * Test max DCS and CF ratios are set correctly
 *
 */
BOOST_AUTO_TEST_CASE(test_max_dcs_ratio)
{
    // These params should give us a,b,c = 9, 252, 45
    FitterUtil::DecayParams_t DecayParams =
        FitterUtil::DecayParams_t{.x = 1, .y = 2, .r = 3, .z_im = 4, .z_re = 5, .width = 6};
    double maxTime = 1;

    // No efficiency
    double efficiencyTimescale = 0;
    auto   cfRate              = [&](double x) { return Phys::cfRate(x, DecayParams, efficiencyTimescale); };
    auto   dcsRate             = [&](double x) { return Phys::dcsRate(x, DecayParams, efficiencyTimescale); };
    auto   gen                 = [](void) { return 1.; };
    auto   genPDF              = [&](double x) {
        return std::exp(-DecayParams.width * x) * DecayParams.width / (1 - std::exp(-DecayParams.width * maxTime));
    };
    std::shared_ptr<std::mt19937> generator;

    SimulatedDecays MyDecays = SimulatedDecays(gen, genPDF, cfRate, dcsRate, std::make_pair(0., maxTime), generator);

    // Apparently the algorithm is so good that float equality works here
    std::cout << MyDecays.maxCFRatio() << std::endl;

    // expect this to be ~0.16, so check it is correct to the machine double epsilon
    BOOST_CHECK(std::abs(MyDecays.maxCFRatio() -
                         1. / (DecayParams.width / (1 - std::exp(-DecayParams.width * maxTime)))) < DBL_EPSILON);

    // expect this to be ~50.8, so check it is correct to ~100x machine double epsilon
    // that should be good enough..
    BOOST_CHECK(
        std::abs(MyDecays.maxDCSRatio() - 306. / (DecayParams.width / (1 - std::exp(-DecayParams.width * maxTime)))) <
        100 * DBL_EPSILON);
}
