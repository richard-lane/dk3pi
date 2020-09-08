#include <cmath>

#include <boost/test/unit_test.hpp>

#include "findZ.h"

/*
 * Test that the right error gets thrown when the avg strong phase is outside of the range -pi to pi
 */
BOOST_AUTO_TEST_CASE(test_phase_out_of_range)
{
    BOOST_CHECK_NO_THROW(amplitudeFiddleFactor(1000, 10, M_PI_2));

    BOOST_CHECK_THROW(amplitudeFiddleFactor(1000, 10, -M_PI - 1), PhaseOutOfRangeException);
    BOOST_CHECK_THROW(amplitudeFiddleFactor(1000, 10, M_PI + 1), PhaseOutOfRangeException);
}

/*
 * Test that the right value is returned by amplitudeFiddleFactor
 */
BOOST_AUTO_TEST_CASE(test_amplitude_fiddle_factor)
{
    size_t numDcs      = 1000;
    size_t numCf       = 10;
    double strongPhase = 1;

    double expectedReal = 10 * std::cos(1);
    double expectedImag = 10 * std::sin(1);

    std::complex<double> fiddleFactor = amplitudeFiddleFactor(numDcs, numCf, strongPhase);
    BOOST_CHECK_CLOSE(expectedReal, fiddleFactor.real(), 1e-4);
    BOOST_CHECK_CLOSE(expectedImag, fiddleFactor.imag(), 1e-4);
}
