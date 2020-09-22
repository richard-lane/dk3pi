#include <boost/test/unit_test.hpp>

#include <cmath>
#include <utility>
#include <vector>

#include "D2K3PiError.h"
#include "RatioCalculator.h"

#define TOLERANCE (1e-8)

/*
 * Check that the ratio of a simple dataset is taken correctly
 */
BOOST_AUTO_TEST_CASE(test_ratio_calculation, *boost::unit_test::tolerance(TOLERANCE))
{
    const std::vector<size_t> myCfDecays{100};
    const std::vector<size_t> myDcsDecays{100};

    // Expect to end up with a vector of ratios (1) and error sqrt(0.02).
    const std::vector<double> expectedRatio{1};
    const std::vector<double> expectedError{std::sqrt(0.02)};

    auto ratioAndError = RatioCalculator::ratioAndError(myCfDecays, myDcsDecays);

    BOOST_CHECK_CLOSE(ratioAndError.first[0], expectedRatio[0], TOLERANCE);
    BOOST_CHECK_CLOSE(ratioAndError.second[0], expectedError[0], TOLERANCE);
}

/*
 * Check an error gets thrown if the datasets are different sizes
 */
BOOST_AUTO_TEST_CASE(test_throw_incompatible_data)
{
    BOOST_CHECK_NO_THROW(RatioCalculator::ratioAndError({1, 2}, {1, 2}));
    BOOST_CHECK_THROW(RatioCalculator::ratioAndError({1, 2, 3}, {1, 2}), RatioCalculator::WrongSizeDataset);
}

/*
 * Check that an error gets thrown if a bin contains zeros
 */
BOOST_AUTO_TEST_CASE(test_throw_zero_points)
{
    const std::vector<size_t> containsZero{1, 0, 1, 1};
    const std::vector<size_t> noZero{1, 1, 1, 1};

    BOOST_CHECK_THROW(RatioCalculator::ratioAndError(containsZero, noZero), RatioCalculator::ZerosInData);
    BOOST_CHECK_THROW(RatioCalculator::ratioAndError(noZero, containsZero), RatioCalculator::ZerosInData);
    BOOST_CHECK_THROW(RatioCalculator::ratioAndError(containsZero, containsZero), RatioCalculator::ZerosInData);
    BOOST_CHECK_NO_THROW(RatioCalculator::ratioAndError(noZero, noZero));
}
