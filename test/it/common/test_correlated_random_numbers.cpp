#include <boost/test/unit_test.hpp>

#include <memory>
#include <random>
#include <vector>

#include "util.h"

BOOST_AUTO_TEST_SUITE(correlated_random_numbers)

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

/*
 * Generate a large-ish set of random numbers, work out the mean/variance/covariance
 */
BOOST_AUTO_TEST_CASE(random_numbers)
{
    // Create a set of random numbers
    size_t                           count = 1e6;
    std::vector<double>              means{3, 5.4, 4.6};
    std::vector<std::vector<double>> cov{std::vector<double>{2, 3.4, -2.8},
                                         std::vector<double>{3.4, 5.84, -4.64},
                                         std::vector<double>{-2.8, -4.64, 4.64}};
    // std::vector<std::vector<double>> cov{
    //    std::vector<double>{1, 0., 0}, std::vector<double>{0., 2, 0}, std::vector<double>{0, 0, 3}};

    std::random_device            rd;
    std::shared_ptr<std::mt19937> gen = std::make_shared<std::mt19937>(rd());

    std::vector<std::vector<double>> correlatedRandomNumbers = util::correlatedGaussianNumbers(gen, count, means, cov);

    // Work out their means, check if they are sane
    std::pair<double, double> stats0 = util::meanAndStdDev(correlatedRandomNumbers[0]);
    std::pair<double, double> stats1 = util::meanAndStdDev(correlatedRandomNumbers[1]);
    std::pair<double, double> stats2 = util::meanAndStdDev(correlatedRandomNumbers[2]);
    std::vector<double>       calculatedMeans{stats0.first, stats1.first, stats2.first};

    for (size_t i = 0; i < 3; i++) {
        BOOST_CHECK(std::fabs(calculatedMeans[i] - means[i]) < 0.01);
    }

    // Work out the covariance matrix, compare it to the other one
    std::vector<std::vector<double>> generatedDataCov = util::covarianceMatrix(correlatedRandomNumbers);
    for (size_t i = 0; i < means.size(); ++i) {
        // 0.5% accuracy is good enough... but it could fail intermittently...
        CHECK_CLOSE_COLLECTIONS(generatedDataCov[i], cov[i], 0.5);
    }
}

BOOST_AUTO_TEST_SUITE_END()
