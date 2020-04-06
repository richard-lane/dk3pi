/*
 * Tests for ratio calculation
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <utility>
#include <vector>

#include "D2K3PiError.h"
#include "RatioCalculator.h"

#define TOLERANCE (0.00000000001)

/*
 * At the moment can't use boose test to check vectors of floats are equal within tolerance;
 * Use this as a workaround
 */
// Have to make it a macro so that it reports exact line numbers when checks fail.
#ifndef CHECK_CLOSE_COLLECTIONS
#define CHECK_CLOSE_COLLECTIONS(aa, bb, tolerance)            \
    {                                                         \
        using std::distance;                                  \
        using std::begin;                                     \
        using std::end;                                       \
        auto a = begin(aa), ae = end(aa);                     \
        auto b = begin(bb);                                   \
        BOOST_CHECK(distance(a, ae) == distance(b, end(bb))); \
        for (; a != ae; ++a, ++b) {                           \
            BOOST_CHECK_CLOSE(*a, *b, tolerance);             \
        }                                                     \
    }
#endif // CHECK_CLOSE_COLLECTIONS

/*
 * Fixture for creating a RatioCalculator
 */
struct TestRatioCalculator {
    TestRatioCalculator(const std::vector<double> &cfDecayTimes,
                        const std::vector<double> &dcsDecayTimes,
                        const std::vector<double> &binLimits)
    {
        Calculator  = new RatioCalculator(cfDecayTimes, dcsDecayTimes, binLimits);
        *Calculator = RatioCalculator(cfDecayTimes, dcsDecayTimes, binLimits);
    }
    ~TestRatioCalculator() { delete Calculator; }
    RatioCalculator *Calculator;
};

/*
 * Check that an exception is raised when unsorted bin limits are passed to the RatioCalculator constructor.
 */
BOOST_AUTO_TEST_CASE(test_unsorted_bin_limits_exception)
{
    const std::vector<double> unsortedVector = {0, 1, 2, 3, 4, 2};
    const std::vector<double> sortedVector   = {1, 2, 3};

    // Create a
    BOOST_CHECK_NO_THROW(TestRatioCalculator bar(unsortedVector, unsortedVector, sortedVector));
    BOOST_CHECK_THROW(TestRatioCalculator foo(sortedVector, sortedVector, unsortedVector), D2K3PiException);
}

/*
 * Check that the ratio of two vectors and its error is calculated correctly
 * Will be comparing floating points for equality, so set a reasonably small value for floating point tolerance.
 */
BOOST_AUTO_TEST_CASE(test_vector_ratio, *boost::unit_test::tolerance(TOLERANCE))
{
    const std::vector<double> vector{0, 1};
    TestRatioCalculator       MyRatioCalculator(vector, vector, vector);

    const std::vector<size_t> numerator{100};
    const std::vector<size_t> denominator{50};

    // Expect error in each to be sqrt(N); combine errors as err = sqrt((N1+N2)/(N1+N2))*ratio
    const double expectedRatio = 2;
    const double expectedError = std::sqrt(0.03) * expectedRatio;

    std::vector<std::pair<double, double>> actualRatioAndError =
        MyRatioCalculator.Calculator->findRatioAndError(numerator, denominator);

    BOOST_CHECK(actualRatioAndError[0].first == expectedRatio);
    BOOST_CHECK(actualRatioAndError[0].second == expectedError);
}

/*
 * Check that the ratio of a simple dataset is taken correctly and that the right attributes are populated.
 */
BOOST_AUTO_TEST_CASE(test_ratio_calculation, *boost::unit_test::tolerance(TOLERANCE))
{
    const std::vector<double> myCfDecays(100, 5);
    const std::vector<double> myDcsDecays(100, 5);
    const std::vector<double> myBinLimits{4, 6};
    TestRatioCalculator       MyRatioCalculator(myCfDecays, myDcsDecays, myBinLimits);

    // Expect to end up with a vector of ratios (1) and error sqrt(0.02).
    const std::vector<double> expectedRatio{1};
    const std::vector<double> expectedError{std::sqrt(0.02)};

    MyRatioCalculator.Calculator->calculateRatios();

    BOOST_CHECK(MyRatioCalculator.Calculator->ratio == expectedRatio);
    BOOST_CHECK(MyRatioCalculator.Calculator->error == expectedError);
}

/*
 * Check that bin centres and widths are set correctly by the constructor
 */
BOOST_AUTO_TEST_CASE(test_bin_centres_widths, *boost::unit_test::tolerance(TOLERANCE))
{
    const std::vector<double> binLimits{-1, 1, 5, 10, 16};
    TestRatioCalculator       MyRatioCalculator(binLimits, binLimits, binLimits);

    // Expect bin centres and widths (0, 3, 7.5, 13) and (2, 4, 5, 6)
    const std::vector<double> expectedBinCentres{0, 3, 7.5, 13};
    const std::vector<double> expectedBinWidths{2, 4, 5, 6};

    BOOST_CHECK(MyRatioCalculator.Calculator->binCentres == expectedBinCentres);
    BOOST_CHECK(MyRatioCalculator.Calculator->binWidths == expectedBinWidths);
}

/*
 * Check that an error gets thrown if a bin contains zeros
 */
BOOST_AUTO_TEST_CASE(test_throw_zero_points)
{
    const std::vector<double> binLimits{0, 1, 2, 3, 4};
    const std::vector<double> containsZero{0.5, 2.5, 3.5};
    const std::vector<double> noZero{0.5, 1.5, 2.5, 3.5};

    TestRatioCalculator MyRatioCalculator1(containsZero, noZero, binLimits);
    TestRatioCalculator MyRatioCalculator2(noZero, containsZero, binLimits);
    TestRatioCalculator MyRatioCalculator3(containsZero, containsZero, binLimits);
    TestRatioCalculator MyRatioCalculator4(noZero, noZero, binLimits);

    BOOST_CHECK_THROW(MyRatioCalculator1.Calculator->calculateRatios(), D2K3PiException);
    BOOST_CHECK_THROW(MyRatioCalculator2.Calculator->calculateRatios(), D2K3PiException);
    BOOST_CHECK_THROW(MyRatioCalculator3.Calculator->calculateRatios(), D2K3PiException);
    BOOST_CHECK_NO_THROW(MyRatioCalculator4.Calculator->calculateRatios());
}

/*
 * Integration-style test that the ratio calculator actually works
 */
BOOST_AUTO_TEST_CASE(it_ratio_calculator, *boost::unit_test::tolerance(TOLERANCE))
{
    // Create a dataset that looks like:
    //     100 50 30 20 15
    //      50 30 20 15  3
    //   points in each bin, then calculate the ratio
    std::vector<double> numerator{};
    std::vector<double> denominator{};
    std::vector<double> binLimits{0, 1, 2, 3, 4, 5};

    // yep
    for (size_t i = 0; i < 100; ++i) {
        numerator.push_back(0.5);
    }

    for (size_t i = 0; i < 50; ++i) {
        numerator.push_back(1.5);
        denominator.push_back(0.6);
    }

    for (size_t i = 0; i < 30; ++i) {
        numerator.push_back(2.4);
        denominator.push_back(1.5);
    }

    for (size_t i = 0; i < 20; ++i) {
        numerator.push_back(3.7);
        denominator.push_back(2.3);
    }

    for (size_t i = 0; i < 15; ++i) {
        numerator.push_back(4.5);
        denominator.push_back(3.6);
    }

    for (size_t i = 0; i < 3; ++i) {
        denominator.push_back(4.7);
    }

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
    std::cout << "\t" << expectedErrors[1] << " " << MyRatioCalculator.error[1] << std::endl;

    CHECK_CLOSE_COLLECTIONS(expectedRatios, MyRatioCalculator.ratio, TOLERANCE);
    CHECK_CLOSE_COLLECTIONS(expectedErrors, MyRatioCalculator.error, TOLERANCE);
}
