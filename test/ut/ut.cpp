// For some reason the ordering here matters
// I think I've accidentally declared the BOOST unit test main() fcn somewhere in one of these

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

// I know
#include "test_util.cpp"
#include "test_ratio_calculator.cpp"
#include "test_fitter.cpp"
#include "test_pull_study.cpp"
#include "test_simulator.cpp"
#include "test_fitter_utils.cpp"
#include "efficiency/test_amplitude.cpp"
#include "efficiency/test_readroot.cpp"
#include "efficiency/test_z.cpp"
#include "efficiency/test_efficiency_util.cpp"
#include "efficiency/test_flat_phsp.cpp"
