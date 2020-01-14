#define BOOST_TEST_MODULE TestHello
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/filesystem/path.hpp>

#include "../include/util.h"
#include "../src/util.cpp"

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
