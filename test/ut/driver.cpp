/*
 * This #define is necessary to define main() for Boost UT
 * 
 * This file should be linked exactly once with each UT executable
 */
#define BOOST_TEST_MODULE "UT"
#include <boost/test/included/unit_test.hpp>
