#ifndef TEST_EFFICIENCY_CPP
#define TEST_EFFICIENCY_CPP

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cstring>

#include <TH1D.h>
#include <TH2D.h>

#include "efficiency.h"
#include "util.h"

/*
 * Test 1d hist entropy
 */
BOOST_AUTO_TEST_CASE(test_entropy)
{
    std::unique_ptr<TH1D> hist = std::make_unique<TH1D>("tmp", "tmp", 5, 0., 5.);
    hist->SetBinContent(1, 0);
    hist->SetBinContent(2, 1);
    hist->SetBinContent(3, 2);
    hist->SetBinContent(4, 2);
    hist->SetBinContent(5, 5);

    double expectedEntropy = 1.220607265;

    BOOST_CHECK_CLOSE(entropy(hist.get()), expectedEntropy, 1e-7);
}

/*
 * Test we can retrieve the right 2d histogram
 */
BOOST_AUTO_TEST_CASE(test_2d_hist)
{
    std::vector<double> bin{0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
    PhspBins            bins{bin, bin, bin, bin, bin};

    PhspPoint              point1{1, 2, 3, 4, 5};
    PhspPoint              point2{5, 4, 3, 2, 1};
    std::vector<PhspPoint> points{point1, point2};

    HistogramProjections binning(bins, points, "test");

    // Check our histograms are the same by comparing their names
    // C-style str comparison => they're the same if strcmp returns 0
    BOOST_CHECK(strcmp(binning.get2dhistogram(1, 2).GetName(), binning.get2dhistogram(2, 1).GetName()) == 0);
}

/*
 * Test the right error is raised if invalid indices are passed to the index finding helper fcn
 */
BOOST_AUTO_TEST_CASE(test_index_conversion_error)
{
    std::vector<double> bin{0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
    PhspBins            bins{bin, bin, bin, bin, bin};

    PhspPoint              point1{1, 2, 3, 4, 5};
    PhspPoint              point2{5, 4, 3, 2, 1};
    std::vector<PhspPoint> points{point1, point2};

    HistogramProjections binning(bins, points, "test");

    BOOST_CHECK_NO_THROW(binning.get1dhistogram(2));
    BOOST_CHECK_THROW(binning.get1dhistogram(5), HistogramNotFound);

    BOOST_CHECK_NO_THROW(binning.get2dhistogram(1, 2));
    BOOST_CHECK_THROW(binning.get2dhistogram(2, 2), HistogramNotFound);
}

/*
 * Test items get binned right by HistogramProjections ctor (and hence the public functions)
 */
BOOST_AUTO_TEST_CASE(test_efficiency_binning)
{
    std::vector<double> bin{0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
    PhspBins            bins{bin, bin, bin, bin, bin, bin};

    PhspPoint              point1{1, 2, 3, 4, 5, 6};
    PhspPoint              point2{5, 4, 3, 2, 1, 6};
    std::vector<PhspPoint> points{point1, point2};

    HistogramProjections binning(bins, points, "test");

    // Expect to have bins like (0, 0, 2, 0, 0) for the p3 bin
    // And (1, 0, 0, 0, 1) for the p0 bin
    TH1D hist2 = binning.get1dhistogram(2);
    // Root bin indexing starts at 1, obviously
    BOOST_CHECK_CLOSE(hist2.GetBinContent(1), 0, 1e-7);
    BOOST_CHECK_CLOSE(hist2.GetBinContent(3), 2, 1e-7);

    TH1D hist0 = binning.get1dhistogram(0);
    BOOST_CHECK_CLOSE(hist0.GetBinContent(2), 0, 1e-7);
    BOOST_CHECK_CLOSE(hist0.GetBinContent(1), 1, 1e-7);

    // Expect to have 1 point in the p0p1 hist at (0, 1), one at (4, 3) and none in (2, 2)
    TH2D hist01 = binning.get2dhistogram(1, 0);
    BOOST_CHECK_CLOSE(hist01.GetBinContent(1, 2), 1, 1e-7);
    BOOST_CHECK_CLOSE(hist01.GetBinContent(5, 4), 1, 1e-7);
    BOOST_CHECK_CLOSE(hist01.GetBinContent(2, 2), 0, 1e-7);
}

/*
 * Test taking the ratio between 2 efficiency binning objects doesn't affect the original objects and returns the right
 * ratio
 */
BOOST_AUTO_TEST_CASE(test_efficiency_binning_ratio)
{
    std::vector<double> bin{0.5, 1.5, 2.5, 3.5};
    PhspBins            bins{bin, bin, bin};

    PhspPoint              point1{1, 1, 1};
    PhspPoint              point2{3, 3, 3};
    std::vector<PhspPoint> hist1Pts{point1, point1, point1, point2};
    std::vector<PhspPoint> hist2Pts{point1, point1, point2, point2, point2};

    HistogramProjections numerator(bins, hist1Pts, "num");
    HistogramProjections denominator(bins, hist2Pts, "denom");

    HistogramProjections ratio = numerator / denominator;

    // Expect ratio of points in 1st 1d hist to be (1.5, 0/0, 1/3)
    TH1D ratioHist1d = ratio.get1dhistogram(0);
    BOOST_CHECK_CLOSE(ratioHist1d.GetBinContent(1), 1.5, 1e-7);
    BOOST_CHECK_CLOSE(ratioHist1d.GetBinContent(2), 0, 1e-7);
    BOOST_CHECK_CLOSE(ratioHist1d.GetBinContent(3), 1 / 3., 1e-7);

    // Expect ratio of points in 2nd 1d hist to be:
    std::array<std::array<double, 3>, 3> expected    = {{{1.5, 0, 0}, {0., 0., 0.}, {0, 0, 1. / 3}}};
    TH2D                                 ratioHist2d = ratio.get2dhistogram(1, 0);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            BOOST_CHECK_CLOSE(ratioHist2d.GetBinContent(i + 1, j + 1), expected[i][j], 1e-7);
        }
    }

    // Check that the numerator and denominator are unaffected
    TH1D numHist = numerator.get1dhistogram(0);
    BOOST_CHECK_CLOSE(numHist.GetBinContent(1), 3, 1e-7);
    BOOST_CHECK_CLOSE(numHist.GetBinContent(2), 0, 1e-7);
    BOOST_CHECK_CLOSE(numHist.GetBinContent(3), 1, 1e-7);

    TH1D denomHist = denominator.get1dhistogram(0);
    BOOST_CHECK_CLOSE(denomHist.GetBinContent(1), 2, 1e-7);
    BOOST_CHECK_CLOSE(denomHist.GetBinContent(2), 0, 1e-7);
    BOOST_CHECK_CLOSE(denomHist.GetBinContent(3), 3, 1e-7);

    TH2D numerator2dHist = numerator.get2dhistogram(0, 1);
    expected             = {{{3, 0, 0}, {0., 0., 0.}, {0, 0, 1}}};
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            BOOST_CHECK_CLOSE(numerator2dHist.GetBinContent(i + 1, j + 1), expected[i][j], 1e-7);
        }
    }
}

/*
 * Test mutual info calculation
 */
BOOST_AUTO_TEST_CASE(test_mutual_info)
{
    // Create a 2D histogram of a dataset with a known, nonzero mutual information
    std::vector<std::vector<size_t>> values = {std::vector<size_t>{4, 2, 1, 1},
                                               std::vector<size_t>{2, 3, 1, 1},
                                               std::vector<size_t>{2, 3, 2, 2},
                                               std::vector<size_t>{8, 0, 0, 0}};
    std::unique_ptr<TH2D>            hist2d(new TH2D("test_mutual_info", "test", 4, 0, 5, 4, 0, 5));
    for (size_t i = 1; i <= values.size(); ++i) {
        for (size_t j = 1; j <= values[i - 1].size(); ++j) {
            hist2d->SetBinContent(i, j, values[i - 1][j - 1]);
        }
    }
    BOOST_CHECK_CLOSE(mutual_info(hist2d.get()), 0.32111859787111086, 1e-7);

    // Create a hist with perfect mutual information
    std::vector<std::vector<size_t>> values_i1 = {std::vector<size_t>{10, 0, 0, 0},
                                                  std::vector<size_t>{0, 8, 0, 0},
                                                  std::vector<size_t>{0, 0, 10, 0},
                                                  std::vector<size_t>{0, 0, 0, 10}};
    for (size_t i = 1; i <= values_i1.size(); ++i) {
        for (size_t j = 1; j <= values_i1[i - 1].size(); ++j) {
            hist2d->SetBinContent(i, j, values_i1[i - 1][j - 1]);
        }
    }
    BOOST_CHECK_CLOSE(mutual_info(hist2d.get()), 1, 1e-7);

    // Create a hist with 0 mutual information
    std::vector<std::vector<size_t>> values_i0 = {std::vector<size_t>{1, 1, 1, 1},
                                                  std::vector<size_t>{1, 1, 1, 1},
                                                  std::vector<size_t>{1, 1, 1, 1},
                                                  std::vector<size_t>{1, 1, 1, 1}};
    for (size_t i = 1; i <= values_i0.size(); ++i) {
        for (size_t j = 1; j <= values_i0[i - 1].size(); ++j) {
            hist2d->SetBinContent(i, j, values_i0[i - 1][j - 1]);
        }
    }
    BOOST_CHECK_CLOSE(mutual_info(hist2d.get()), 0, 1e-7);
}

/*
 * Test swapping histogram axes
 */
BOOST_AUTO_TEST_CASE(test_swap_2d_hist)
{
    TH2D originalHist = TH2D("title", "name", 2, 0, 1, 2, 0, 1);
    originalHist.SetBinContent(1, 1, 1);
    originalHist.SetBinContent(1, 2, 2);
    originalHist.SetBinContent(2, 1, 3);
    originalHist.SetBinContent(2, 2, 4);

    TH2D swappedHist = swapAxes(originalHist);

    BOOST_CHECK_CLOSE(swappedHist.GetBinContent(1, 1), 1, 1e-7);
    BOOST_CHECK_CLOSE(swappedHist.GetBinContent(1, 2), 3, 1e-7);
    BOOST_CHECK_CLOSE(swappedHist.GetBinContent(2, 1), 2, 1e-7);
    BOOST_CHECK_CLOSE(swappedHist.GetBinContent(2, 2), 4, 1e-7);
}

#endif // TEST_EFFICIENCY_CPP
