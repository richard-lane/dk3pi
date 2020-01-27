#ifndef TEST_FITTER_CPP
#define TEST_FITTER_CPP

#define BOOST_TEST_DYN_LINK
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <stdio.h>
#include <vector>

#include "../../include/D2K3PiError.h"
#include "../../include/Fitter.h"
#include "../../src/Fitter.cpp"

/*
 * At the moment can't use boose test to check vectors of floats are equal within tolerance;
 * Use this as a workaround
 */
// Have to make it a macro so that it reports exact line numbers when checks fail.
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

/*
 * Test that passing in data and bins of differing lengths causes an error
 */
BOOST_AUTO_TEST_CASE(test_bad_fitdata)
{
    std::vector<double> threeLength = {1, 2, 3};
    std::vector<double> twoLength   = {1, 2};

    BOOST_CHECK_THROW(FitData MyData(threeLength, threeLength, twoLength, twoLength), D2K3PiException);
}

/*
 * Test that creating a FitData object with overlapping bin limits causes an error
 */
BOOST_AUTO_TEST_CASE(test_overlapping_bins)
{
    std::vector<double> linspace{1, 2, 3};
    std::vector<double> widths{1.5, 1.5, 1.5};

    BOOST_CHECK_THROW(FitData MyData(linspace, widths, linspace, linspace), D2K3PiException);
}

/*
 * Test that creating a FitData object containing inf causes an error
 */
BOOST_AUTO_TEST_CASE(test_inf_data_error)
{
    std::vector<double> v      = {1, 2};
    std::vector<double> widths = {0.1, 0.1};
    std::vector<double> infs   = {1.0 / 0.0, 1.0 / 0.0};

    BOOST_CHECK_THROW(FitData MyData(v, widths, infs, v), D2K3PiException);
}

/*
 * Test that creating a FitData object containing NaNs causes an error
 */
BOOST_AUTO_TEST_CASE(test_nan_data_error)
{
    std::vector<double> v      = {1, 2};
    std::vector<double> widths = {0.1, 0.1};
    std::vector<double> nans   = {0.0 / 0.0, 0.0 / 0.0};

    BOOST_CHECK_THROW(FitData MyData(v, widths, nans, v), D2K3PiException);
}

/*
 * Test that creating a FitData object containing zeros causes an error
 */
BOOST_AUTO_TEST_CASE(test_zero_data_error)
{
    std::vector<double> v      = {1, 2};
    std::vector<double> widths = {0.1, 0.1};
    std::vector<double> zeros  = {0.0, 0.0};

    BOOST_CHECK_THROW(FitData MyData(v, widths, zeros, v), D2K3PiException);
}

/*
 * Test that unsorted bin centres causes an error
 */
BOOST_AUTO_TEST_CASE(test_unsorted_bin_centres_exception)
{
    std::vector<double> centres{1, 2, 1.5};
    std::vector<double> widths{0.1, 0.1, 0.1};

    BOOST_CHECK_THROW(FitData MyData(centres, widths, widths, widths), D2K3PiException);
}

/*
 * Test that data following a quadratic function has sensible fit parameters
 *
 * Again this is more like IT than UT but no one is going to notice
 */
BOOST_AUTO_TEST_CASE(test_quadratic_fit, *boost::unit_test::tolerance(0.01))
{
    // Set our x data to linear spaced bins of width 0.1 centred on integers
    std::vector<double> binCentres{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<double> binWidths(11, 0.1);

    // Use the function 1 + x + x^2
    // Set the errors to something reasonable but small
    std::vector<double> data{1, 3, 7, 13, 21, 31, 43, 57, 73, 91, 111};
    std::vector<double> errors(11, 0.01);

    // We expect the fit coefficients to all be 1
    std::vector<double> expectedFitCoefficients{1, 1, 1};

    // Create a fitter with our parameters
    FitData_t fitData(binCentres, binWidths, data, errors);
    Fitter    MyFitter(fitData);

    // Perform a fit and check that our coefficients are all 1
    MyFitter.pol2fit();
    CHECK_CLOSE_COLLECTIONS(MyFitter.fitParams.fitParams, expectedFitCoefficients, 0.001);
}

/*
 * Test that a plot cannot be made without running the fitter first.
 */
BOOST_AUTO_TEST_CASE(test_plot_before_fit)
{
    std::vector<double> oneTwo{1, 2};
    std::vector<double> zeros{0, 0};
    FitData_t           fitData(oneTwo, zeros, oneTwo, oneTwo);
    Fitter              MyFitter(fitData);

    BOOST_CHECK_THROW(MyFitter.saveFitPlot("foo", "bar"), D2K3PiException);
}

/*
 * Test that we cannot write a plot to a file that already exists
 */
BOOST_AUTO_TEST_CASE(test_plot_already_exists)
{
    std::vector<double> oneTwo{1, 2};
    std::vector<double> zeros{0, 0};
    FitData_t           fitData(oneTwo, zeros, oneTwo, oneTwo);
    Fitter              MyFitter(fitData);

    // Change the fit parameter to a non empty vector so that we will not hit the plot before fit error
    MyFitter.fitParams.fitParams = {1};

    // Create a temporary file and attempt to save a plot there
    // Cleanup isn't automatic, so this might leave temp files all over the place
    boost::filesystem::path temp   = boost::filesystem::unique_path();
    std::string             native = temp.native();
    std::ofstream(temp.native());
    BOOST_CHECK_THROW(MyFitter.saveFitPlot("foo", native), D2K3PiException);

    // Remove the temporary file
    remove(native.c_str());
}

/*
 * Test that a file is created when saveFitPlot is called
 * Doesn't check that the plot is sensible- just checks that it is created.
 */
BOOST_AUTO_TEST_CASE(test_plot_created)
{
    std::vector<double> oneTwo{1, 2};
    std::vector<double> zeros{0, 0};
    FitData_t           fitData(oneTwo, zeros, oneTwo, oneTwo);
    Fitter              MyFitter(fitData);

    // Change the fit parameter to a non empty vector so that we will not hit the plot before fit error
    MyFitter.fitParams.fitParams = {1};

    // Create a temporary file and attempt to save a plot there
    boost::filesystem::path temp   = boost::filesystem::unique_path();
    std::string             native = temp.native();
    MyFitter.saveFitPlot("foo", native);

    // Remove the temporary file
    remove(native.c_str());
}

#endif // TEST_FITTER_CPP
