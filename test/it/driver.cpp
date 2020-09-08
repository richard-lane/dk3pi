/*
 * This #define is necessary to define main() for Boost testing
 * 
 * This file should be linked exactly once with each IT executable
 */
#define BOOST_TEST_MODULE "IT"
#include <boost/test/included/unit_test.hpp>
