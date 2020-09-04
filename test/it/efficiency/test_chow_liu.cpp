#define BOOST_TEST_DYN_LINK
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <TH3D.h>

#include "ChowLiu.h"
#include "util.h"

BOOST_AUTO_TEST_SUITE(test_chow_liu)

/*
 * Create a 3d histogram with a probability of 0.1 to find an event in a given bin and see if we get something close
 * back
 */
BOOST_AUTO_TEST_CASE(test_chow_liu_probability)
{
    PhspBins               bins{{0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}};
    ChowLiu::Approximation Distribution(bins, "test");

    // Put 1 point in each 3d bin
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 3; ++k) {
                PhspPoint point{i + 0.5, j + 0.5, k + 0.5};
                Distribution.binPoint(point);
            }
        }
    }

    Distribution.makeApproximation();
    BOOST_CHECK_CLOSE(Distribution.value({1.5, 0.5, 0.5}), 1 / 27., 1e-7);
}

/*
 * Create 3d Gaussian, check the approximation gives us the right thing
 */
BOOST_AUTO_TEST_CASE(test_chow_liu_gaussians)
{
    const size_t dimension{3};
    const size_t numBins{100};
    const double min{-4.};
    const double max{4.};

    // Define uniformly spaced bins
    PhspBins bins(dimension);
    for (auto& bin : bins) {
        bin = std::vector<double>(numBins + 1);
        for (size_t i = 0; i <= numBins; ++i) {
            bin[i] = min + i * (max - min) / numBins;
        }
    }

    // Create ROOT histogram and an approximation object thing
    TH3D rootHist("root Hist", "root Hist", numBins, bins[0].data(), numBins, bins[1].data(), numBins, bins[2].data());
    ChowLiu::Approximation Approximation(bins, "approximation");

    // Generate points from a 3d gaussian and populate the ROOT histogram and the approximation object thing
    const size_t                     numDataPoints{1000000};
    std::random_device               rd;
    std::shared_ptr<std::mt19937>    gen = std::make_shared<std::mt19937>(rd());
    std::vector<double>              means{0, 1, -1};
    std::vector<std::vector<double>> covarianceMatrix{{2, 0.5, 0.3}, {0.5, 0.5, -0.2}, {0.3, -0.2, 1}};
    std::vector<std::vector<double>> randomNumbers =
        util::correlatedGaussianNumbers(gen, numDataPoints, means, covarianceMatrix);
    for (size_t i = 0; i < numDataPoints; ++i) {
        rootHist.Fill(randomNumbers[0][i], randomNumbers[1][i], randomNumbers[2][i]);
        Approximation.binPoint({randomNumbers[0][i], randomNumbers[1][i], randomNumbers[2][i]});
    }

    // Create a hist for reconstructed bin contents
    TH3D approxHist("approx", "approx", numBins, bins[0].data(), numBins, bins[1].data(), numBins, bins[2].data());
    Approximation.makeApproximation();

    // Iterate over every bin in our approximate histogram, finding the Chow-Liu approximate probability in each one
    // Then scale that up to give the number of points in each bin predicted by our approximation
    for (size_t i = 0; i < numBins; ++i) {
        double x = (bins[0][i] + bins[0][i + 1]) / 2;

        for (size_t j = 0; j < numBins; ++j) {
            double y = (bins[1][j] + bins[1][j + 1]) / 2;

            for (size_t k = 0; k < numBins; ++k) {
                double z = (bins[2][k] + bins[2][k + 1]) / 2;

                // Find the chow liu estimate for probability and work out how many events we'd expect in that bin
                double prob = Approximation.value({x, y, z});
                approxHist.SetBinContent(rootHist.FindBin(x, y, z), prob * numDataPoints);
            }
        }
    }
    util::saveObjectToFile(&rootHist, "root.png");
    util::saveObjectToFile(&approxHist, "approx.png");

    // Check that the probability that our histograms are the same is 1 (ish)
    // i.e. check that our histograms are the same
    double pValue = approxHist.Chi2Test(&rootHist, "P");
    BOOST_CHECK_CLOSE(pValue, 1, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
