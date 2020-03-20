#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include "MinuitFitter.h"

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

    double expectedChiSq = 5.395436507936507;

    // Create a PolynomialFitFcn object and find chi squared
    SimplePolynomialFunction MyPolynomial(1.0, 1.0, 1.0);
    BasePolynomialFcn        MyFitFcn(data, times, errors);

    BOOST_CHECK(MyFitFcn.operator()(std::vector<double>{1.0, 1.0, 1.0}) == expectedChiSq);
}
