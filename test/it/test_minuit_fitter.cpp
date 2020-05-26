#include <boost/test/unit_test.hpp>

#include <vector>

#include "FitterUtils.h"
#include "MinuitPolynomialFitter.h"
#include "PhysicalFitter.h"

BOOST_AUTO_TEST_SUITE(fitter_IT)

/*
 * Test chi squared comes out right when finding chi sq, not integrating
 */
BOOST_AUTO_TEST_CASE(test_chi_squared, *boost::unit_test::tolerance(0.000000000000001))
{
    // Create vectors representing our data and errors
    // Theoretically, this dataset might fit to 1 + t + t^2
    std::vector<double> times{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    std::vector<double> data{2, 10, 15, 20, 31, 40, 55, 70, 90, 113};
    std::vector<double> errors{2, 4, 6, 8, 10, 12, 14, 16, 18, 20};

    double expectedChiSq = 1.070386944129503;

    // Create a PolynomialFitFcn object and find chi squared
    IntegralOptions_t  integralOptions(false, 0, errors, 0);
    PolynomialChiSqFcn MyFitFcn(data, times, errors, integralOptions);

    BOOST_CHECK(std::fabs(MyFitFcn(std::vector<double>{1.0, 1.0, 1.0}) - expectedChiSq) < 1e-7);
}

/*
 * Test constraints are implemented properly
 */
BOOST_AUTO_TEST_CASE(test_constraints)
{
    // Create a fitter object such that chi squared just evaluates to the constraint
    std::vector<double> data{};
    std::vector<double> binLimits{0};
    IntegralOptions_t   integralOptions(true, 0, binLimits, 0);
    FitData_t           fitData(binLimits, data, data);
    PhysicalFitter      fitter(fitData, integralOptions, true);

    // Fix all parameters except X and Y
    fitter.setFitParams(std::vector<double>{1, 1, 1, 1, 1, 1}, std::vector<double>{1, 1, 1, 1, 1, 1});
    fitter.fixParameters(std::vector<std::string>{"z_re", "z_im", "width", "r"});

    // Minimise
    fitter.fit();

    // Check that we have recovered the constraint...
    BOOST_CHECK_CLOSE(fitter.fitParams.fitParams[0], WORLD_AVERAGE_X, 1e-7);
    BOOST_CHECK_CLOSE(fitter.fitParams.fitParams[1], WORLD_AVERAGE_Y, 1e-7);
    BOOST_CHECK_CLOSE(fitter.fitParams.correlationMatrix->GetMatrixArray()[1], X_Y_CORRELATION, 1e-7);
}

BOOST_AUTO_TEST_SUITE_END()
