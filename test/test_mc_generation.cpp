/*
 * Tests for MC generation
 */
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/filesystem/path.hpp>
#include <utility>

#include "../include/D2K3PiError.h"
#include "../include/MCGenerator.h"
#include "../src/MCGenerator.cpp"

// ---- Fixtures
/*
 * Monte Carlo generator produced with the default range
 */
struct DefaultGenerator {
    DefaultGenerator()
    {
        Generator  = new MCGenerator;
        *Generator = MCGenerator();
    }
    ~DefaultGenerator() { delete Generator; }
    MCGenerator *Generator;
};

/*
 * Monte Carlo generator produced with a range between 5 and 15.
 */
struct SpecialGenerator {
    SpecialGenerator(double a, double b, double c, double d)
    {
        Generator  = new MCGenerator;
        *Generator = MCGenerator(std::make_pair(a, b), std::make_pair(c, d));
    }
    ~SpecialGenerator() { delete Generator; }
    MCGenerator *Generator;
};

/*
 * A quadratic function y = x^2 + 1
 */
const double quadratic(double x)
{
    return (x * x) + 1;
}

// ---- Test cases
// TODO these distribution tests should really be IT
/*
 * Check that the default x values follow the expected distribution.
 */
BOOST_AUTO_TEST_CASE(test_default_generated_x_vals)
{
    DefaultGenerator Gen;
    double           tot = 0;
    for (int i = 0; i < 10000; ++i) {
        tot += Gen.Generator->getRandomX();
    }

    // Expect the average randomly generated number to be 0.5
    // This range is huge!
    BOOST_CHECK(tot < 5500 && tot > 4500);
}

/*
 * Check that x values are between the expected range
 */
BOOST_AUTO_TEST_CASE(test_specified_generated_x_vals)
{
    SpecialGenerator Gen(5, 15, 0, 1);
    double           tot = 0;
    for (int i = 0; i < 10000; ++i) {
        tot += Gen.Generator->getRandomX();
    }

    // Expect the average randomly generated number to be 10
    // This range is huge so won't catch all errors :(
    BOOST_CHECK(tot < 105000 && tot > 95000);
}

/*
 * Check that the default y values follow the expected distribution.
 */
BOOST_AUTO_TEST_CASE(test_default_generated_y_vals)
{
    DefaultGenerator Gen;
    double           tot = 0;
    for (int i = 0; i < 10000; ++i) {
        tot += Gen.Generator->getRandomY();
    }

    // Expect the average randomly generated number to be 0.5
    // This range is huge!
    BOOST_CHECK(tot < 5500 && tot > 4500);
}

/*
 * Check that y values are between the expected range
 */
BOOST_AUTO_TEST_CASE(test_specified_generated_y_vals)
{
    SpecialGenerator Gen(0, 1, 3, 17);
    double           tot = 0;
    for (int i = 0; i < 10000; ++i) {
        tot += Gen.Generator->getRandomY();
    }

    // Expect the average randomly generated number to be 10
    // This range is huge so won't catch all errors :(
    BOOST_CHECK(tot < 105000 && tot > 95000);
}

/*
 * Check that (x, y) is rejected from the distribution f(x) with y > f(x)
 */
BOOST_AUTO_TEST_CASE(test_rejection)
{
    SpecialGenerator Gen(0, 1, 0, 2);
    BOOST_CHECK(!Gen.Generator->isAccepted(0.5, 1.5, &quadratic));
}

/*
 * Check that (x, y) is accepted as being from the distribution f(x) if y < f(x)
 */
BOOST_AUTO_TEST_CASE(test_acceptance)
{
    SpecialGenerator Gen(0, 1, 0, 2);
    BOOST_CHECK(Gen.Generator->isAccepted(1, 1, &quadratic));
}

/*
 * Check that an appropriate exception is raised when the generated function evaluates to a value outside of the allowed
 * region.
 */
BOOST_AUTO_TEST_CASE(test_function_out_of_range_exception)
{
    DefaultGenerator Gen;
    BOOST_CHECK_THROW(Gen.Generator->isAccepted(0.75, 1, &quadratic), D2K3PiException);
}

/*
 * Check that an appropriate exception is raised when the provided x value is outside of the allowed range.
 */
BOOST_AUTO_TEST_CASE(test_x_out_of_range_exception)
{
    SpecialGenerator Gen(0, 1, 0, 5);

    // Must be careful that our implementation catches the right error while we don't have specific error classes
    // Provided x value is smaller than the allowed minimum
    BOOST_CHECK_THROW(Gen.Generator->isAccepted(-0.5, 1, &quadratic), D2K3PiException);

    // Provided x value is larger than the allowed maximum
    BOOST_CHECK_THROW(Gen.Generator->isAccepted(1.5, 1, &quadratic), D2K3PiException);
}

/*
 * Check that an appropriate exception is raised when the provided y value is outside of the allowed range.
 */
BOOST_AUTO_TEST_CASE(test_y_out_of_range_exception)
{
    SpecialGenerator Gen(0, 1, 0, 5);

    // Must be careful that our implementation catches the right error while we don't have specific error classes
    // Provided x value is smaller than the allowed minimum
    BOOST_CHECK_THROW(Gen.Generator->isAccepted(0.5, -1, &quadratic), D2K3PiException);

    // Provided x value is larger than the allowed maximum
    BOOST_CHECK_THROW(Gen.Generator->isAccepted(0.5, 6, &quadratic), D2K3PiException);
}

/*
 * Check that a malformed std::pair is detected
 */
BOOST_AUTO_TEST_CASE(test_x_pair_malformed_exception)
{
    SpecialGenerator Gen(0, 1, 0, 5);

    BOOST_CHECK_THROW(Gen.Generator->setXRange(std::make_pair(1, 0)), D2K3PiException);
}

BOOST_AUTO_TEST_CASE(test_y_pair_malformed_exception)
{
    SpecialGenerator Gen(0, 1, 0, 5);

    BOOST_CHECK_THROW(Gen.Generator->setYRange(std::make_pair(1, 0)), D2K3PiException);
}