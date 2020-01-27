/*
 * Tests for ratio calculation
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <utility>
#include <vector>

#include "../../include/D2K3PiError.h"
#include "../../include/RatioCalculator.h"
#include "../../src/RatioCalculator.cpp"

#define TOLERANCE (0.00000000001)

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
 * Check that a vector is binned as expected.
 */
BOOST_AUTO_TEST_CASE(test_vector_binning)
{
    const std::vector<double> decayTimes = {1, 6, 3, 4, 5, 2, 7, 8, 9};
    const std::vector<double> binLimits  = {0.5, 1.5, 3.5, 7.5, 9.5};

    TestRatioCalculator MyRatioCalculator(decayTimes, decayTimes, binLimits);

    // With this data expect bins of size (1, 2, 4, 2)
    std::vector<size_t> expectedNumPerBin = {1, 2, 4, 2};

    BOOST_CHECK(MyRatioCalculator.Calculator->binVector(decayTimes, binLimits) == expectedNumPerBin);
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
 * Check that bad ratios are pruned from the data, i.e. that ratios of 0/0, 1/0 or 0/1 are removed but 1/1 remain.
 * Also check that the bins corresponding to these are removed from binCentres and binWidths.
 *
 * Not really a unit test, because this uses the whole workflow and doesn't isolate the bit of code that we're
 * interested in, but we can just pretend that this was meant to be IT and it's somehow better to have more IT and less
 * UT
 */
BOOST_AUTO_TEST_CASE(test_prune_data, *boost::unit_test::tolerance(TOLERANCE))
{
    // This data should give us 4 bins with (0, 1, 0, 1) CF points and (1, 0, 0, 1) DCS points.
    // We can thus check that only the final bin remains after pruning.
    const std::vector<double> cfTimes{1.5, 3.5};
    const std::vector<double> dcsTimes{0.5, 3.6};
    const std::vector<double> binLimits{0, 1, 2, 3, 4};

    TestRatioCalculator       MyRatioCalculator(cfTimes, dcsTimes, binLimits);
    const std::vector<double> expectedRatios{1};
    const std::vector<double> expectedErrors{std::sqrt(2)};
    const std::vector<double> expectedBinCentres{3.5};
    const std::vector<double> expectedBinWidths{1};

    // When calculating ratios, we expect that the bad ratios are pruned.
    // This isn't really testing in isolation => isn't really a unit test, but it's good enough...
    MyRatioCalculator.Calculator->calculateRatios();
    BOOST_CHECK(MyRatioCalculator.Calculator->binCentres == expectedBinCentres);
    BOOST_CHECK(MyRatioCalculator.Calculator->binWidths == expectedBinWidths);
    BOOST_CHECK(MyRatioCalculator.Calculator->ratio == expectedRatios);
    BOOST_CHECK(MyRatioCalculator.Calculator->error == expectedErrors);
}
