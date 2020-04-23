#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include "MinuitFcns.h"

/*
 * Test chi squared comes out right
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
    PolynomialChiSqFcnNoIntegral MyFitFcn(data, times, errors);

    BOOST_CHECK(std::fabs(MyFitFcn(std::vector<double>{1.0, 1.0, 1.0}) - expectedChiSq) < 1e-7);
}
