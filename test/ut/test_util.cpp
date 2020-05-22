#define BOOST_TEST_MODULE UT
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/filesystem/path.hpp>
#include <cfloat>

#include "D2K3PiError.h"
#include "physics.h"
#include "util.h"

#include "TGraph.h"

BOOST_AUTO_TEST_CASE(test_find_bin_limits_unsorted_data)
{
    std::vector<double> unsortedData = {1, 2, 3, 4, 2.5};

    BOOST_CHECK_THROW(util::findBinLimits(unsortedData, 1, 0, 5), D2K3PiException);
}

BOOST_AUTO_TEST_CASE(test_find_bin_limits_too_few_points)
{
    std::vector<double> data = {1, 2, 3, 4, 5};

    BOOST_CHECK_THROW(util::findBinLimits(data, 6, 0, 6), D2K3PiException);
}

BOOST_AUTO_TEST_CASE(test_find_bin_limits_low_edge_too_high)
{
    std::vector<double> data = {1, 2, 3, 4, 5};
    BOOST_CHECK_THROW(util::findBinLimits(data, 2, 1.5, 6), D2K3PiException);
}

BOOST_AUTO_TEST_CASE(test_find_bin_limits_high_edge_too_low)
{
    std::vector<double> data = {1, 2, 3, 4, 5};
    BOOST_CHECK_THROW(util::findBinLimits(data, 2, 0, 4.5), D2K3PiException);
}

BOOST_AUTO_TEST_CASE(test_0_maximum_points)
{
    std::vector<double> data = {1, 2, 3};
    BOOST_CHECK_THROW(util::findBinLimits(data, 0, 0.5, 3.5), D2K3PiException);
}

BOOST_AUTO_TEST_CASE(test_find_bin_limits, *boost::unit_test::tolerance(10 * DBL_EPSILON))
{
    std::vector<double> data              = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<double> expectedBinLimits = {0.5, 3.5, 6.5, 10.5};
    std::vector<double> actualBinLimits   = util::findBinLimits(data, 3, 0.5, 10.5);
    BOOST_CHECK(expectedBinLimits == actualBinLimits);
}

/*
 * Check that a vector is binned as expected.
 */
BOOST_AUTO_TEST_CASE(test_vector_binning)
{
    const std::vector<double> decayTimes = {1, 6, 3, 4, 5, 2, 7, 8, 9};
    const std::vector<double> binLimits  = {0.5, 1.5, 3.5, 7.5, 9.5};

    // With this data expect bins of size (1, 2, 4, 2)
    std::vector<size_t> expectedNumPerBin = {1, 2, 4, 2};

    BOOST_CHECK(util::binVector(decayTimes, binLimits) == expectedNumPerBin);
}

/*
 * Check that passing vectors with different numbers of objects/options or 0 options causes an err
 */
BOOST_AUTO_TEST_CASE(test_draw_multiple_objects)
{
    TGraph *MyTGraph1 = new TGraph();
    TGraph *MyTGraph2 = new TGraph();

    const std::vector<TObject *> twoTGraphs{MyTGraph1, MyTGraph2};
    const std::vector<TObject *> oneTGraph{MyTGraph1};
    const std::vector<TObject *> zeroTGraphs{};

    const std::vector<std::string> oneString{"a"};
    const std::vector<std::string> zeroStrings{};

    const util::LegendParams_t legendParams;
    std::vector<std::string>   legendLabels{"a", "b"};

    BOOST_CHECK_THROW(util::saveObjectsToFile<TGraph>(twoTGraphs, oneString, legendLabels, "path", legendParams),
                      D2K3PiException);
    BOOST_CHECK_THROW(util::saveObjectsToFile<TGraph>(zeroTGraphs, oneString, legendLabels, "path", legendParams),
                      D2K3PiException);
    BOOST_CHECK_THROW(util::saveObjectsToFile<TGraph>(twoTGraphs, zeroStrings, legendLabels, "path", legendParams),
                      D2K3PiException);
    BOOST_CHECK_THROW(util::saveObjectsToFile<TGraph>(zeroTGraphs, zeroStrings, legendLabels, "path", legendParams),
                      D2K3PiException);

    delete MyTGraph1;
    delete MyTGraph2;
}

/*
 * Test rate integral calculators
 */
BOOST_AUTO_TEST_CASE(test_integral_calculators, *boost::unit_test::tolerance(1e-8))
{
    DecayParams_t DecayParams = {
        .x = 0.0, .y = std::sqrt(0.12), .r = 1, .z_im = 0.0, .z_re = 0.2 / std::sqrt(0.12), .width = 10};

    // Cursory check that expectedParams still works
    BOOST_CHECK(std::abs(util::expectedParams(DecayParams)[0] - 1) < 1e-8);
    BOOST_CHECK(std::abs(util::expectedParams(DecayParams)[1] - 2) < 1e-8);
    BOOST_CHECK(std::abs(util::expectedParams(DecayParams)[2] - 3) < 1e-8);

    // DCS integrals
    // Use efficiency timescale of 0 such that the efficiency is unity
    BOOST_CHECK(std::abs(Phys::dcsIntegralWithEfficiency(
                             0, 3, util::expectedParams(DecayParams), DecayParams.width, 0, 1e-14, 25) -
                         0.12599999999966256411) < 1e-10);

    // CF integrals
    // Use efficiency timescale of 0 such that the efficiency is unity
    BOOST_CHECK(std::abs(Phys::cfIntegralWithEfficiency(0, 3, DecayParams.width, 0, 1e-15, 25) -
                         0.09999999999999064237703) < 1e-10);
}

/*
 * Test decay rates
 */
BOOST_AUTO_TEST_CASE(test_rates, *boost::unit_test::tolerance(1e-8))
{
    DecayParams_t DecayParams = {.x = 1, .y = 2, .r = 3, .z_im = 4, .z_re = 5, .width = 6};

    // Trust that util::expectedParams works
    std::vector<double> expectedParams = util::expectedParams(DecayParams);

    // Rate ratio
    BOOST_CHECK_SMALL(Phys::rateRatio(100, expectedParams) - 475209.0, 1e-10);
    BOOST_CHECK_SMALL(Phys::rateRatio(100, DecayParams) - 475209.0, 1e-10);

    // DCS rate
    BOOST_CHECK_SMALL(Phys::wrongSignDecayRate(2, DecayParams) - 693 * std::exp(-DecayParams.width * 2), 1e-10);
    BOOST_CHECK_SMALL(
        Phys::wrongSignDecayRate(2, expectedParams, DecayParams.width) - 693 * std::exp(-DecayParams.width * 2), 1e-10);

    // CF rate
    BOOST_CHECK_SMALL(Phys::rightSignDecayRate(2, DecayParams) - std::exp(-DecayParams.width * 2), 1e-10);
    BOOST_CHECK_SMALL(Phys::rightSignDecayRate(2, DecayParams.width) - std::exp(-DecayParams.width * 2), 1e-10);
}
