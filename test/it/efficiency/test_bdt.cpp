#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <random>

#include "bdt_reweighting.h"
#include "efficiencyUtil.h"
#include "util.h"

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL BDT_NP_ARRAY_API
#include <numpy/arrayobject.h>

BOOST_AUTO_TEST_SUITE(test_bdt)

/*
 * Test that a simple efficiency (0.5 everywhere) gets recovered
 */
BOOST_AUTO_TEST_CASE(test_bdt_simple_efficiency)
{
    // mcData is a 10^5 point Normal distribution
    // realData is the same distribution with half the number of points
    size_t                 numPoints = 100000;
    std::vector<PhspPoint> mcData(numPoints, std::vector<double>(1));
    std::vector<PhspPoint> realData(numPoints / 2, std::vector<double>(1));

    std::random_device               rd;
    std::mt19937                     gen(rd());
    std::normal_distribution<double> gaussian;
    for (auto& point : mcData) {
        // Each point is a 1-element vector
        point[0] = gaussian(gen);
    }
    for (auto& point : mcData) {
        point[0] = gaussian(gen);
    }

    // Pass these data to the BDT to get an estimate of our efficiency function
    PyObject* bdt = initBDT(mcData, realData);

    // Create another Gaussian and find the weights associated with each point
    std::vector<PhspPoint> moreData(realData.size(), std::vector<double>(1));
    for (auto& point : moreData) {
        point[0] = gaussian(gen);
    }
    std::vector<double> bdtEfficiencies{efficiency(bdt, moreData, numPoints)};

    // The average weight should be almost exactly 2
    BOOST_CHECK_CLOSE(util::findMean(bdtEfficiencies), 2.0, 0.01);
}

BOOST_AUTO_TEST_SUITE_END()
