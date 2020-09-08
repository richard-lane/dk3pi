#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include "D2K3PiError.h"
#include "FitterUtils.h"

/*
 * Test that passing in data and bins of differing lengths causes an error
 */
BOOST_AUTO_TEST_CASE(test_bad_fitdata)
{
    std::vector<double> twoLength = {1, 2};

    BOOST_CHECK_THROW(FitData MyData(twoLength, twoLength, twoLength), D2K3PiException);
}

/*
 * Test that unsorted bin limits cause an error
 */
BOOST_AUTO_TEST_CASE(test_unsorted_bin_centres_exception)
{
    std::vector<double> limits{1, 2, 1.5};
    std::vector<double> data{0., 0.};
    BOOST_CHECK_THROW(FitData MyData(limits, data, data), D2K3PiException);
}

/*
 * Test that creating a FitData object with nearly overlapping bins doesn't cause an error
 */
BOOST_AUTO_TEST_CASE(test_nearly_overlapping_bins)
{
    std::vector<double> timeBinLimits{0, 1, 2, 3, 4};
    std::vector<double> data{1, 1, 1, 1};

    // Pass our stuff into FitData; the widths will be close but shouldn't cause an error.
    BOOST_CHECK_NO_THROW(FitData MyData(timeBinLimits, data, data));
}

/*
 * Test that creating a FitData object containing inf causes an error
 */
BOOST_AUTO_TEST_CASE(test_inf_data_error)
{
    std::vector<double> binLimits = {0, 1, 2};
    std::vector<double> infs      = {1.0 / 0.0, 1.0 / 0.0};
    std::vector<double> errs      = {0.0, 0.0};

    BOOST_CHECK_THROW(FitData MyData(binLimits, infs, errs), D2K3PiException);
}

/*
 * Test that creating a FitData object containing NaNs causes an error
 */
BOOST_AUTO_TEST_CASE(test_nan_data_error)
{
    std::vector<double> nans      = {0.0 / 0.0, 0.0 / 0.0};
    std::vector<double> binLimits = {0, 1, 2};
    std::vector<double> errs      = {0.0, 0.0};

    BOOST_CHECK_THROW(FitData MyData(binLimits, nans, errs), D2K3PiException);
}

/*
 * Test that creating a FitData object containing zeros causes an error
 */
BOOST_AUTO_TEST_CASE(test_zero_data_error)
{
    std::vector<double> zeros     = {0.0, 0.0};
    std::vector<double> binLimits = {0, 1, 2};
    std::vector<double> errs      = {0.0, 0.0};

    BOOST_CHECK_THROW(FitData MyData(binLimits, zeros, errs), D2K3PiException);
}
