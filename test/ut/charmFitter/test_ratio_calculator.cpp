#include <boost/test/unit_test.hpp>

#include <cmath>
#include <utility>
#include <vector>

#include "D2K3PiError.h"
#include "RatioCalculator.h"

#define TOLERANCE (0.00000000001)

/*
 * Check that an exception is raised when unsorted bin limits are passed to the RatioCalculator constructor.
 */
BOOST_AUTO_TEST_CASE(test_unsorted_bin_limits_exception)
{
    const std::vector<double> unsortedVector = {0, 1, 2, 3, 4, 2};
    const std::vector<double> sortedVector   = {1, 2, 3};
    const std::vector<size_t> counts{};

    BOOST_CHECK_NO_THROW(RatioCalculator bar(counts, counts, sortedVector));
    BOOST_CHECK_THROW(RatioCalculator foo(counts, counts, unsortedVector), D2K3PiException);
}

/*
 * Check that the ratio of a simple dataset is taken correctly and that the right attributes are populated.
 */
BOOST_AUTO_TEST_CASE(test_ratio_calculation, *boost::unit_test::tolerance(TOLERANCE))
{
    const std::vector<size_t> myCfDecays{100};
    const std::vector<size_t> myDcsDecays{100};
    const std::vector<double> myBinLimits{4, 6};
    RatioCalculator           MyRatioCalculator(myCfDecays, myDcsDecays, myBinLimits);

    // Expect to end up with a vector of ratios (1) and error sqrt(0.02).
    const std::vector<double> expectedRatio{1};
    const std::vector<double> expectedError{std::sqrt(0.02)};

    MyRatioCalculator.calculateRatios();

    BOOST_CHECK(MyRatioCalculator.ratio == expectedRatio);
    BOOST_CHECK(MyRatioCalculator.error == expectedError);
}

/*
 * Check that bin centres and widths are set correctly by the constructor
 */
BOOST_AUTO_TEST_CASE(test_bin_centres_widths, *boost::unit_test::tolerance(TOLERANCE))
{
    const std::vector<double> binLimits{-1, 1, 5, 10, 16};
    const std::vector<size_t> counts{};
    RatioCalculator           MyRatioCalculator(counts, counts, binLimits);

    // Expect bin centres and widths (0, 3, 7.5, 13) and (2, 4, 5, 6)
    const std::vector<double> expectedBinCentres{0, 3, 7.5, 13};
    const std::vector<double> expectedBinWidths{2, 4, 5, 6};

    BOOST_CHECK(MyRatioCalculator.binCentres == expectedBinCentres);
    BOOST_CHECK(MyRatioCalculator.binWidths == expectedBinWidths);
}

/*
 * Check that an error gets thrown if a bin contains zeros
 */
BOOST_AUTO_TEST_CASE(test_throw_zero_points)
{
    const std::vector<double> binLimits{0, 1, 2, 3, 4};
    const std::vector<size_t> containsZero{1, 0, 1, 1};
    const std::vector<size_t> noZero{1, 1, 1, 1};

    RatioCalculator MyRatioCalculator1(containsZero, noZero, binLimits);
    RatioCalculator MyRatioCalculator2(noZero, containsZero, binLimits);
    RatioCalculator MyRatioCalculator3(containsZero, containsZero, binLimits);
    RatioCalculator MyRatioCalculator4(noZero, noZero, binLimits);

    BOOST_CHECK_THROW(MyRatioCalculator1.calculateRatios(), D2K3PiException);
    BOOST_CHECK_THROW(MyRatioCalculator2.calculateRatios(), D2K3PiException);
    BOOST_CHECK_THROW(MyRatioCalculator3.calculateRatios(), D2K3PiException);
    BOOST_CHECK_NO_THROW(MyRatioCalculator4.calculateRatios());
}
