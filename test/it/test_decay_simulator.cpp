#ifndef TEST_DECAY_SIMULATOR_CPP
#define TEST_DECAY_SIMULATOR_CPP

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>

BOOST_AUTO_TEST_CASE(none)
{
    std::cout << "none" << std::endl;
    BOOST_CHECK(true);
}

#endif // TEST_DECAY_SIMULATOR_CPP
