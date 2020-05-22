#ifndef TEST_RATIO_CALCULATOR
#define TEST_RATIO_CALCULATOR

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include "RatioCalculator.h"

#include <cmath>
/*
 * Integration-style test that the ratio calculator actually works
 */
BOOST_AUTO_TEST_CASE(it_ratio_calculator)
{
    // Create a dataset and calculate the ratio
    std::vector<size_t> numerator{100, 50, 30, 20, 15};
    std::vector<size_t> denominator{50, 30, 20, 15, 3};
    std::vector<double> binLimits{0, 1, 2, 3, 4, 5};

    // Find our expected ratios and their errors
    std::vector<double> expectedRatios{2.0, 5.0 / 3.0, 3.0 / 2.0, 4.0 / 3.0, 5.0};
    std::vector<double> expectedErrors{expectedRatios[0] * std::sqrt(0.03),
                                       expectedRatios[1] * std::sqrt(4.0 / 75.0),
                                       expectedRatios[2] * std::sqrt(1.0 / 12.0),
                                       expectedRatios[3] * std::sqrt(7.0 / 60.0),
                                       expectedRatios[4] * std::sqrt(2.0 / 5.0)};

    // Use the ratio calculator to bin times and calculate ratios
    RatioCalculator MyRatioCalculator = RatioCalculator(denominator, numerator, binLimits);
    MyRatioCalculator.calculateRatios();

    BOOST_CHECK_CLOSE(expectedRatios[0], MyRatioCalculator.ratio[0], 1e-10);
    BOOST_CHECK_CLOSE(expectedRatios[1], MyRatioCalculator.ratio[1], 1e-10);
    BOOST_CHECK_CLOSE(expectedRatios[2], MyRatioCalculator.ratio[2], 1e-10);
    BOOST_CHECK_CLOSE(expectedRatios[3], MyRatioCalculator.ratio[3], 1e-10);
    BOOST_CHECK_CLOSE(expectedRatios[4], MyRatioCalculator.ratio[4], 1e-10);

    BOOST_CHECK_CLOSE(expectedErrors[0], MyRatioCalculator.error[0], 1e-10);
    BOOST_CHECK_CLOSE(expectedErrors[1], MyRatioCalculator.error[1], 1e-10);
    BOOST_CHECK_CLOSE(expectedErrors[2], MyRatioCalculator.error[2], 1e-10);
    BOOST_CHECK_CLOSE(expectedErrors[3], MyRatioCalculator.error[3], 1e-10);
    BOOST_CHECK_CLOSE(expectedErrors[4], MyRatioCalculator.error[4], 1e-10);
}

#endif // TEST_RATIO_CALCULATOR
