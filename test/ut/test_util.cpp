#define BOOST_TEST_MODULE UT
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/filesystem/path.hpp>
#include <cfloat>

#include "../../include/D2K3PiError.h"
#include "../../include/util.h"
#include "../../src/util.cpp"

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
