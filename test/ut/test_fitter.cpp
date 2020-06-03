#ifndef TEST_FITTER_CPP
#define TEST_FITTER_CPP

#define BOOST_TEST_DYN_LINK
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <stdio.h>
#include <vector>

#include "D2K3PiError.h"
#include "RatioCalculator.h"
#include "fitter/FitterUtils.h"
#include "fitter/MinuitPolynomialFitter.h"
#include "fitter/PhysicalFitter.h"
#include "fitter/RootFitter.h"

#include "TMatrixD.h"

/*
 * Test that data following a quadratic function has sensible fit parameters
 *
 * Again this is more like IT than UT but no one is going to notice
 */
BOOST_AUTO_TEST_CASE(test_quadratic_fit, *boost::unit_test::tolerance(0.01))
{
    // Set our x data to linear spaced bins centred on integers
    std::vector<double> binLimits{-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5};

    // Use the function 1 + x + x^2
    // Set the errors to something reasonable but small
    std::vector<double> data{1, 3, 7, 13, 21, 31, 43, 57, 73, 91, 111};
    std::vector<double> errors(11, 0.01);

    // We expect the fit coefficients to all be 1
    std::vector<double> expectedFitCoefficients{1, 1, 1};

    // Create a fitter with our parameters
    FitData_t  fitData(binLimits, data, errors);
    RootFitter MyFitter(fitData);

    // Perform a fit and check that our coefficients are all 1
    MyFitter.fit(0, 1, "Q");
    CHECK_CLOSE_COLLECTIONS(MyFitter.fitParams.fitParams, expectedFitCoefficients, 0.001);
}

/*
 * Test that data following a quadratic function has sensible fit parameters, when we fit using Minuit2 chisq
 *
 * this is more like IT than UT but no one is going to notice
 */
BOOST_AUTO_TEST_CASE(test_quadratic_fit_minuit2_chi_sq, *boost::unit_test::tolerance(0.0001))
{
    // Set our x data to linear spaced bins centred on integers
    std::vector<double> binLimits{-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5};

    // Use the function 1 + x + x^2
    // Set the errors to something reasonable but small
    std::vector<double> data{1, 3, 7, 13, 21, 31, 43, 57, 73, 91, 111};
    std::vector<double> errors(11, 0.01);

    // We expect the fit coefficients to all be 1
    std::vector<double> expectedFitCoefficients{1, 1, 1};

    // Create a fitter with our parameters
    FitData_t              fitData(binLimits, data, errors);
    IntegralOptions_t      integralOptions(false, 0, errors, 0);
    MinuitPolynomialFitter MyFitter(fitData, integralOptions);

    // Perform a fit and check that our coefficients are all 1
    std::vector<double> parameterGuess{0.9, 1.1, 1.2};
    std::vector<double> errorGuess{1, 1, 1};
    MyFitter.setPolynomialParams(parameterGuess, errorGuess);
    MyFitter.fit();
    CHECK_CLOSE_COLLECTIONS(MyFitter.fitParams.fitParams, expectedFitCoefficients, 0.001);
}

/*
 * Test that a plot cannot be made without running the fitter first.
 */
BOOST_AUTO_TEST_CASE(test_plot_before_fit)
{
    std::vector<double> oneTwoThree{1, 2, 3};
    std::vector<double> oneTwo{1, 2};
    FitData_t           fitData(oneTwoThree, oneTwo, oneTwo);
    RootFitter          MyFitter(fitData);

    BOOST_CHECK_THROW(MyFitter.saveFitPlot("foo", "bar"), D2K3PiException);
}

/*
 * Test that we cannot write a plot to a file that already exists
 */
BOOST_AUTO_TEST_CASE(test_plot_already_exists)
{
    std::vector<double> oneTwoThree{1, 2, 3};
    std::vector<double> oneTwo{1, 2};
    FitData_t           fitData(oneTwoThree, oneTwo, oneTwo);
    RootFitter          MyFitter(fitData);

    // Run the fit, which allows us to attempt to save the plot
    MyFitter.fit(0, 1, "Q");

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
 * This definitely isn't UT
 */
BOOST_AUTO_TEST_CASE(test_plot_created)
{
    std::vector<double> oneTwoThree{1, 2, 3};
    std::vector<double> oneTwo{1, 2};
    FitData_t           fitData(oneTwoThree, oneTwo, oneTwo);
    RootFitter          MyFitter(fitData);

    // Run the fit, which allows us to attempt to save the plot
    MyFitter.fit(0, 1, "Q");

    // Create a temporary file and attempt to save a plot there
    boost::filesystem::path temp   = boost::filesystem::unique_path();
    std::string             native = temp.native();
    MyFitter.saveFitPlot("foo", native);
    BOOST_CHECK(boost::filesystem::exists(native));

    // Remove the temporary file
    remove(native.c_str());
}

/*
 * Test that the correlation matrix gets assigned correctly
 */
BOOST_AUTO_TEST_CASE(test_matrix_assigned)
{
    std::vector<double> oneTwoThree{1, 2, 3};
    std::vector<double> oneTwo{1, 2};
    FitData_t           fitData(oneTwoThree, oneTwo, oneTwo);
    RootFitter          MyFitter(fitData);

    // Run the fit, which allows us to attempt to save the plot
    MyFitter.fit(0, 1, "Q");

    // Check that the correlation matrix has been assigned some memory
    BOOST_CHECK(MyFitter.fitParams.correlationMatrix != nullptr);

    // Check that the correlation matrix can be used
    BOOST_CHECK_NO_THROW(MyFitter.fitParams.correlationMatrix->Print());
}

/*
 * Test that our conversion from a vector representing a Covariance matrix to a correlation TMatrixD works
 */
BOOST_AUTO_TEST_CASE(test_corr_cov_conversion, *boost::unit_test::tolerance(0.000000001))
{
    // The correlation matrix we expect
    // Minuit's covariance matrix is represented as a vector (a b c d e f ...):
    //     a
    //     b c
    //     d e f
    // etc.
    std::vector<double> covarianceVector{4, 2, 100, 3, 5, 9};
    std::vector<double> errors{2, 10, 3};
    std::vector<double> expectedCorrMatrixValues{1, 0.1, 0.5, 0.1, 1, 1. / 6., 0.5, 1. / 6., 1};
    TMatrixD            expectedCorrMatrix = TMatrixD(3, 3, expectedCorrMatrixValues.data());

    // Create a fitter object and set the errors to the right things. We will need to pass sensible bin limits to the
    // fitter
    std::vector<double>    binLimits{0, 1, 2, 3};
    FitData_t              fitData(binLimits, errors, errors);
    IntegralOptions_t      integralOptions(false, 0, errors, 0);
    MinuitPolynomialFitter MyFitter(fitData, integralOptions);
    MyFitter.setPolynomialParams(std::vector<double>(3, 1), errors);

    // Hack here where I set something that shouldn't be possible
    MyFitter.fitParams.fitParamErrors = errors;

    BOOST_CHECK(MyFitter.covarianceVector2CorrelationMatrix(covarianceVector) == expectedCorrMatrix);
}

#endif // TEST_FITTER_CPP
