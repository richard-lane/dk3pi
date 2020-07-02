#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <TFile.h>

#include "ReadRoot.h"
#include "amplitudes.h"

BOOST_AUTO_TEST_SUITE(test_amplitude)

/*
 * Test that the correct error gets thrown when we can't find the specified fcn in a shared library
 *
 * NB this will fail if the test is not run in <d2k3pi dir>/build/ or a similar path, which I cba to fix
 */
BOOST_AUTO_TEST_CASE(test_amplitude)
{
    // Read root files to find decay parameters
    boost::filesystem::path currentDir = boost::filesystem::path(__FILE__).parent_path();
    boost::filesystem::path cfPath("cf.root");
    TFile*                  cfFile = new TFile((currentDir / cfPath).string().c_str());

    boost::filesystem::path dcsPath("dcs.root");
    TFile*                  dcsFile = new TFile((currentDir / dcsPath).string().c_str());

    // Branch names
    std::vector<std::string> branchNames = {"_1_K~", "_2_pi#", "_3_pi#", "_4_pi~"};

    ReadRoot cf(cfFile, "DalitzEventList", branchNames);
    ReadRoot dcs(dcsFile, "DalitzEventList", branchNames);

    // Find DCS and CF amplitude functions
    boost::filesystem::path cfLibraryPath("../../../AmpGenTools/amplitude_models/cf.so");
    boost::filesystem::path dcsLibraryPath("../../../AmpGenTools/amplitude_models/dcs.so");
    boost::filesystem::path cfLib  = currentDir / cfLibraryPath;
    boost::filesystem::path dcsLib = currentDir / dcsLibraryPath;

    amplitudeFcnPtr dcsFunc = readFromSharedLib(dcsLib.string());
    amplitudeFcnPtr cfFunc  = readFromSharedLib(cfLib.string());

    std::cout << "DCS decays:"
              << "\n";
    std::cout << "DCS component\t\tCF component" << std::endl;
    for (size_t i = 0; i < 10; ++i) {
        // Find DCS/CF amplitudes of decays that were generated to be DCS
        std::complex<double> dcsComponentDcs = amplitude(dcs.events[i], dcsFunc);
        std::complex<double> dcsComponentCf  = amplitude(dcs.events[i], cfFunc);

        // Print components of "DCS" decays
        std::cout << dcsComponentDcs << "\t";
        std::cout << dcsComponentCf << std::endl;
    }

    std::cout << "\nCF decays:"
              << "\n";
    std::cout << "DCS component\t\tCF component" << std::endl;
    for (size_t i = 0; i < 10; ++i) {
        // Find DCS/CF amplitudes of decays that were generated to be CF
        std::complex<double> cfComponentDcs = amplitude(cf.events[i], dcsFunc);
        std::complex<double> cfComponentCf  = amplitude(cf.events[i], cfFunc);

        // Print components of "DCS" decays
        std::cout << cfComponentDcs << "\t";
        std::cout << cfComponentCf << std::endl;
    }

    delete cfFile;
    delete dcsFile;
}

BOOST_AUTO_TEST_SUITE_END()
