#define BOOST_TEST_MODULE UT
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/filesystem/path.hpp>
#include <cfloat>

#include "D2K3PiError.h"
#include "util.h"

#include "TGraph.h"

BOOST_AUTO_TEST_CASE(test_concat_paths_no_extension)
{
    std::string expectedPath{"foo/bar"};
    BOOST_CHECK(util::concatPaths("foo", "bar", "") == expectedPath);
}

BOOST_AUTO_TEST_CASE(test_concat_paths_file_extension)
{
    std::string expectedPath{"foo/bar.pdf"};
    BOOST_CHECK(util::concatPaths("foo", "bar", ".pdf") == expectedPath);
}

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
 * Check that passing vectors with different numbers of objects/options or 0 options causes an err
 */
BOOST_AUTO_TEST_CASE(test_draw_multiple_objects_)
{
    TGraph *MyTGraph1 = new TGraph();
    TGraph *MyTGraph2 = new TGraph();

    const std::vector<TObject *> twoTGraphs{MyTGraph1, MyTGraph2};
    const std::vector<TObject *> oneTGraph{MyTGraph1};
    const std::vector<TObject *> zeroTGraphs{};

    const std::vector<std::string> oneString{"a"};
    const std::vector<std::string> zeroStrings{};

    BOOST_CHECK_THROW(util::saveObjectsToFile<TGraph>(twoTGraphs, oneString, "path"), D2K3PiException);
    BOOST_CHECK_THROW(util::saveObjectsToFile<TGraph>(zeroTGraphs, oneString, "path"), D2K3PiException);
    BOOST_CHECK_THROW(util::saveObjectsToFile<TGraph>(twoTGraphs, zeroStrings, "path"), D2K3PiException);
    BOOST_CHECK_THROW(util::saveObjectsToFile<TGraph>(zeroTGraphs, zeroStrings, "path"), D2K3PiException);

    delete MyTGraph1;
    delete MyTGraph2;
}
