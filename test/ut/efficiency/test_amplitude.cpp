#include <boost/test/unit_test.hpp>

#include "amplitudes.h"

#include <boost/filesystem.hpp>

/*
 * Test that the correct error gets thrown when we can't find the specified fcn in a shared library
 *
 * NB this will fail if the test is not run in <d2k3pi dir>/build/ or a similar path, which I cba to fix
 */
BOOST_AUTO_TEST_CASE(test_amplitude_fcn_not_found)
{
    boost::filesystem::path currentDir = boost::filesystem::path(__FILE__).parent_path();
    boost::filesystem::path libraryPath("../../../AmpGenTools/amplitude_models/cf.so");
    boost::filesystem::path library = currentDir / libraryPath;

    std::string fake_fcn = "NOT_REAL_FUNCTION";
    std::string real_fcn = "AMP";

    BOOST_CHECK_NO_THROW(readFromSharedLib(library.string(), real_fcn));
    BOOST_CHECK_THROW(readFromSharedLib(library.string(), fake_fcn), FunctionNotFoundException);
}
