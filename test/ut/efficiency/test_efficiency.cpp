#ifndef TEST_EFFICIENCY_CPP
#define TEST_EFFICIENCY_CPP

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <TH1D.h>
#include <TH2D.h>

#include "efficiency.h"
#include "util.h"

/*
 * Test mutual info calculation
 */
BOOST_AUTO_TEST_CASE(test_mutual_info)
{
    // Create a 2D histogram of a dataset where the binning variables contain no mutual information

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

#endif // TEST_EFFICIENCY_CPP
