#include <boost/test/unit_test.hpp>

#include <boost/filesystem/path.hpp>
#include <cfloat>

#include "D2K3PiError.h"
#include "fitterUtil.h"
#include "physics.h"
#include "util.h"

#include "TGraph.h"
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

BOOST_AUTO_TEST_CASE(test_find_bin_limits_unsorted_data)
{
    std::vector<double> unsortedData = {1, 2, 3, 4, 2.5};

    BOOST_CHECK_THROW(util::findBinLimits(unsortedData, 1, 0, 5), D2K3PiException);
}

BOOST_AUTO_TEST_CASE(test_find_bin_limits_too_few_points)
{
    std::vector<double> data = {1, 2, 3, 4, 5};

    BOOST_CHECK_THROW(util::findBinLimits(data, 6, 0, 6), D2K3PiException);
}

BOOST_AUTO_TEST_CASE(test_find_bin_limits_low_edge_too_high)
{
    std::vector<double> data = {1, 2, 3, 4, 5};
    BOOST_CHECK_THROW(util::findBinLimits(data, 2, 1.5, 6), D2K3PiException);
}

BOOST_AUTO_TEST_CASE(test_find_bin_limits_high_edge_too_low)
{
    std::vector<double> data = {1, 2, 3, 4, 5};
    BOOST_CHECK_THROW(util::findBinLimits(data, 2, 0, 4.5), D2K3PiException);
}

BOOST_AUTO_TEST_CASE(test_0_maximum_points)
{
    std::vector<double> data = {1, 2, 3};
    BOOST_CHECK_THROW(util::findBinLimits(data, 0, 0.5, 3.5), D2K3PiException);
}

BOOST_AUTO_TEST_CASE(test_find_bin_limits, *boost::unit_test::tolerance(10 * DBL_EPSILON))
{
    std::vector<double> data              = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<double> expectedBinLimits = {0.5, 3.5, 6.5, 10.5};
    std::vector<double> actualBinLimits   = util::findBinLimits(data, 3, 0.5, 10.5);
    BOOST_CHECK(expectedBinLimits == actualBinLimits);
}

/*
 * Check that a vector is binned as expected.
 */
BOOST_AUTO_TEST_CASE(test_vector_binning)
{
    const std::vector<double> decayTimes = {1, 6, 3, 4, 5, 2, 7, 8, 9};
    const std::vector<double> binLimits  = {0.5, 1.5, 3.5, 7.5, 9.5};

    // With this data expect bins of size (1, 2, 4, 2)
    std::vector<size_t> expectedNumPerBin = {1, 2, 4, 2};

    BOOST_CHECK(util::binVector(decayTimes, binLimits) == expectedNumPerBin);
}

/*
 * Test bin overflow
 */
BOOST_AUTO_TEST_CASE(test_bin_overflow)
{
    const std::vector<double> decayTimes = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    const std::vector<double> binLimits  = {0.5, 1.5, 3.5, 7.5, 9.5};

    BOOST_CHECK_THROW(util::binVector(decayTimes, binLimits), D2K3PiException);
}

/*
 * Test bin underflow
 */
BOOST_AUTO_TEST_CASE(test_bin_underflow)
{
    const std::vector<double> decayTimes = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    const std::vector<double> binLimits  = {0.5, 1.5, 3.5, 7.5, 9.5};

    BOOST_CHECK_THROW(util::binVector(decayTimes, binLimits), D2K3PiException);
}

/*
 * Check that passing vectors with different numbers of objects/options or 0 options causes an err
 */
BOOST_AUTO_TEST_CASE(test_draw_multiple_objects)
{
    TGraph *MyTGraph1 = new TGraph();
    TGraph *MyTGraph2 = new TGraph();

    const std::vector<TObject *> twoTGraphs{MyTGraph1, MyTGraph2};
    const std::vector<TObject *> oneTGraph{MyTGraph1};
    const std::vector<TObject *> zeroTGraphs{};

    const std::vector<std::string> oneString{"a"};
    const std::vector<std::string> zeroStrings{};

    const util::LegendParams_t legendParams;
    std::vector<std::string>   legendLabels{"a", "b"};

    BOOST_CHECK_THROW(util::saveObjectsToFile<TGraph>(twoTGraphs, oneString, legendLabels, "path", legendParams),
                      D2K3PiException);
    BOOST_CHECK_THROW(util::saveObjectsToFile<TGraph>(zeroTGraphs, oneString, legendLabels, "path", legendParams),
                      D2K3PiException);
    BOOST_CHECK_THROW(util::saveObjectsToFile<TGraph>(twoTGraphs, zeroStrings, legendLabels, "path", legendParams),
                      D2K3PiException);
    BOOST_CHECK_THROW(util::saveObjectsToFile<TGraph>(zeroTGraphs, zeroStrings, legendLabels, "path", legendParams),
                      D2K3PiException);

    delete MyTGraph1;
    delete MyTGraph2;
}

/*
 * Test we find the correct values of a, b and c
 * Use a small value of tolerance; we are comparing resonably large floats
 * Should move this to test util
 */
BOOST_AUTO_TEST_CASE(test_expected_params, *boost::unit_test::tolerance(0.0000000001))
{
    // Set our decay parameters and what we expect a, b and c to evaluate to
    FitterUtil::DecayParams_t DecayParams =
        FitterUtil::DecayParams_t{.x = 1, .y = 2, .r = 3, .z_im = 4, .z_re = 5, .width = 6};

    // Parameters are described by eq. 2.10 in 1412.7254v2
    double expectedA = 9;
    double expectedB = 252;
    double expectedC = 45;

    auto params = Phys::expectedParams(DecayParams);

    BOOST_CHECK(expectedA == params[0]);
    BOOST_CHECK(expectedB == params[1]);
    BOOST_CHECK(expectedC == params[2]);
}

/*
 * Test rate integrals
 */
BOOST_AUTO_TEST_CASE(test_integrals, *boost::unit_test::tolerance(1e-8))
{
    FitterUtil::DecayParams_t DecayParams = {
        .x = 0.0, .y = std::sqrt(0.12), .r = 1, .z_im = 0.0, .z_re = 0.2 / std::sqrt(0.12), .width = 10};

    // Cursory check that expectedParams still works
    BOOST_CHECK(std::abs(Phys::expectedParams(DecayParams)[0] - 1) < 1e-8);
    BOOST_CHECK(std::abs(Phys::expectedParams(DecayParams)[1] - 2) < 1e-8);
    BOOST_CHECK(std::abs(Phys::expectedParams(DecayParams)[2] - 3) < 1e-8);

    // DCS integrals
    // Use efficiency timescale of 0 such that the efficiency is unity
    BOOST_CHECK(
        std::abs(Phys::dcsIntegralWithEfficiency(0, 3, Phys::expectedParams(DecayParams), DecayParams.width, 0) -
                 0.12599999999966256411) < 1e-10);

    // CF integrals
    // Use efficiency timescale of 0 such that the efficiency is unity
    BOOST_CHECK(std::abs(Phys::cfIntegralWithEfficiency(0, 3, DecayParams.width, 0) - 0.09999999999999064237703) <
                1e-10);
}

/*
 * Test decay rates
 */
BOOST_AUTO_TEST_CASE(test_rates, *boost::unit_test::tolerance(1e-8))
{
    FitterUtil::DecayParams_t DecayParams = {.x = 1, .y = 2, .r = 3, .z_im = 4, .z_re = 5, .width = 6};

    // Trust that Phys::expectedParams works
    auto expectedParams = Phys::expectedParams(DecayParams);

    // Rate ratio
    BOOST_CHECK_SMALL(Phys::rateRatio(100, expectedParams) - 475209.0, 1e-10);
    BOOST_CHECK_SMALL(Phys::rateRatio(100, DecayParams) - 475209.0, 1e-10);

    // DCS rate
    BOOST_CHECK_SMALL(Phys::dcsRate(2, DecayParams, 0) - 693 * std::exp(-DecayParams.width * 2), 1e-10);
    BOOST_CHECK_SMALL(Phys::dcsRate(2, expectedParams, DecayParams.width, 0) - 693 * std::exp(-DecayParams.width * 2),
                      1e-10);

    // CF rate
    BOOST_CHECK_SMALL(Phys::cfRate(2, DecayParams, 0) - std::exp(-DecayParams.width * 2), 1e-10);
    BOOST_CHECK_SMALL(Phys::cfRate(2, DecayParams.width, 0) - std::exp(-DecayParams.width * 2), 1e-10);
}

/*
 * Test conversion between Re/Im and magnitude/phase
 */
BOOST_AUTO_TEST_CASE(test_mag_phase)
{
    std::vector<double> real{2, -1, -1, 2};
    std::vector<double> imaginary{1, 1, -2, -2};

    std::vector<std::pair<double, double>> expectedMagPhase{std::make_pair(2.2361, 0.4636),
                                                            std::make_pair(1.4142, 2.3562),
                                                            std::make_pair(2.2361, 4.2487),
                                                            std::make_pair(2.8284, 5.4978)};
    std::vector<std::pair<double, double>> result = util::reIm2magPhase(real, imaginary);

    BOOST_CHECK_SMALL(result[0].first - expectedMagPhase[0].first, 1e-4);
    BOOST_CHECK_SMALL(result[0].second - expectedMagPhase[0].second, 1e-4);

    BOOST_CHECK_SMALL(result[1].first - expectedMagPhase[1].first, 1e-4);
    BOOST_CHECK_SMALL(result[1].second - expectedMagPhase[1].second, 1e-4);

    BOOST_CHECK_SMALL(result[2].first - expectedMagPhase[2].first, 1e-4);
    BOOST_CHECK_SMALL(result[2].second - expectedMagPhase[2].second, 1e-4);

    BOOST_CHECK_SMALL(result[3].first - expectedMagPhase[3].first, 1e-4);
    BOOST_CHECK_SMALL(result[3].second - expectedMagPhase[3].second, 1e-4);
}

/*
 * Test Gauss-Legendre quadrature with a lambda
 */
BOOST_AUTO_TEST_CASE(test_gauss_legendre_lambda)
{
    auto f = [](double x) { return x * x * std::exp(-x); };

    double low  = 0;
    double high = 4;
    double expectedIntegral =
        (low * low + 2 * low + 2) * std::exp(-low) - (high * high + 2 * high + 2) * std::exp(-high);

    double actualIntegral = FitterUtil::gaussLegendreQuad(f, low, high);

    BOOST_CHECK_SMALL(actualIntegral - expectedIntegral, 1e-15);
}

/*
 * Test covariance matrix
 */
BOOST_AUTO_TEST_CASE(test_cov_matrix)
{
    std::vector<std::vector<double>> data{std::vector<double>{1, 2, 3, 4, 5}, std::vector<double>{2, 4, 5, 7, 9}};

    std::vector<std::vector<double>> expectedCov{std::vector<double>{2, 3.4}, std::vector<double>{3.4, 5.84}};
    std::vector<std::vector<double>> cov = util::covarianceMatrix(data);

    for (size_t i = 0; i < data.size(); ++i) {
        CHECK_CLOSE_COLLECTIONS(expectedCov[i], cov[i], 1e-8);
    }
}

/*
 * Test mean
 */
BOOST_AUTO_TEST_CASE(test_mean)
{
    std::array<double, 10> arr{2, 3, 7, 5, 6, 7, 8, 9, 10, 11};
    BOOST_CHECK_CLOSE(util::findMean(arr), 6.8, 1e-9);
}

/*
 * Test mean + std dev function
 */
BOOST_AUTO_TEST_CASE(test_mean_std)
{
    std::vector<double>       data{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    std::pair<double, double> expected = std::make_pair(6.5, 3.45205);
    std::pair<double, double> actual   = util::meanAndStdDev(data);

    BOOST_CHECK_SMALL(expected.first - actual.first, 1e-5);
    BOOST_CHECK_SMALL(expected.second - actual.second, 1e-5);
}

/*
 * Test matrix multiplication
 */
BOOST_AUTO_TEST_CASE(test_matrix_multiply)
{
    std::vector<std::vector<double>> matrix{std::vector<double>{1, 0, 2, 0},
                                            std::vector<double>{0, 3, 0, 4},
                                            std::vector<double>{0, 0, 5, 0},
                                            std::vector<double>{6, 0, 0, 7}};
    std::vector<double>              vector{2, 5, 1, 8};
    std::vector<double>              expectedResult{4, 47, 5, 68};
    std::vector<double>              result = util::multiply(matrix, vector);

    CHECK_CLOSE_COLLECTIONS(result, expectedResult, 1e-8);
}

/*
 * Test cholesky decomposition with non-square matrix
 */
BOOST_AUTO_TEST_CASE(test_cholesky_decomp_non_square)
{
    std::vector<std::vector<double>> nonSquareMatrix{std::vector<double>{1, 2, 3}, std::vector<double>{4, 5, 6}};
    std::vector<std::vector<double>> malformedMatrix{
        std::vector<double>{1, 2, 3}, std::vector<double>{4, 5, 6}, std::vector<double>{7, 8}};

    BOOST_CHECK_THROW(util::choleskyDecomp(nonSquareMatrix), D2K3PiException);
    BOOST_CHECK_THROW(util::choleskyDecomp(malformedMatrix), D2K3PiException);
}

/*
 * Test cholesky decomposition with non-symmetric matrix
 */
BOOST_AUTO_TEST_CASE(test_cholesky_decomp_non_symmetric)
{
    std::vector<std::vector<double>> asymmetricMatrix{
        std::vector<double>{1, 2, 3}, std::vector<double>{4, 5, 6}, std::vector<double>{7, 8, 9}};

    BOOST_CHECK_THROW(util::choleskyDecomp(asymmetricMatrix), D2K3PiException);
}

/*
 * Test cholesky decomposition
 */
BOOST_AUTO_TEST_CASE(test_cholesky_decomp)
{
    std::vector<std::vector<double>> matrix{
        std::vector<double>{4, 12, -16}, std::vector<double>{12, 37, -43}, std::vector<double>{-16, -43, 98}};
    std::vector<std::vector<double>> expectedLowerTriangular{
        std::vector<double>{2, 0, 0}, std::vector<double>{6, 1, 0}, std::vector<double>{-8, 5, 3}};

    std::vector<std::vector<double>> cholDecomp = util::choleskyDecomp(matrix);

    for (size_t i = 0; i < matrix.size(); ++i) {
        CHECK_CLOSE_COLLECTIONS(cholDecomp[i], expectedLowerTriangular[i], 1e-8);
    }
}

/*
 * Test cholesky decomposition with a 4x4 matrix, just in case
 */
BOOST_AUTO_TEST_CASE(test_cholesky_decomp_4x4)
{
    std::vector<std::vector<double>> matrix{std::vector<double>{18, 22, 54, 42},
                                            std::vector<double>{22, 70, 86, 62},
                                            std::vector<double>{54, 86, 174, 134},
                                            std::vector<double>{42, 62, 134, 106}};
    std::vector<std::vector<double>> expectedLowerTriangular{std::vector<double>{4.24264, 0, 0, 0},
                                                             std::vector<double>{5.18545, 6.56591, 0, 0},
                                                             std::vector<double>{12.72792, 3.04604, 1.64974, 0},
                                                             std::vector<double>{9.89949, 1.62455, 1.84971, 1.39262}};

    std::vector<std::vector<double>> cholDecomp = util::choleskyDecomp(matrix);

    for (size_t i = 0; i < matrix.size(); ++i) {
        CHECK_CLOSE_COLLECTIONS(cholDecomp[i], expectedLowerTriangular[i], 1e-3);
    }
}

/*
 * Test that generating correlated gaussians doesnt throw
 */
BOOST_AUTO_TEST_CASE(test_correlated_gaussian_no_throw)
{
    std::vector<std::vector<double>> covariance{std::vector<double>{1, 0.5}, std::vector<double>{0.5, 1}};
    std::vector<double>              means{1, 2};
    std::shared_ptr<std::mt19937>    generator = std::make_shared<std::mt19937>(0);

    BOOST_CHECK_NO_THROW(util::correlatedGaussianNumbers(generator, 1, means, covariance));
}
