#define BOOST_TEST_DYN_LINK
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <TFile.h>

#include <numeric>

#include "ReadRoot.h"
#include "amplitudes.h"

BOOST_AUTO_TEST_SUITE(test_amplitude)

/*
 * Test that the correct error gets thrown when we can't find the specified fcn in a shared library
 *
 * This doesn't really actually test for much since I don't know what amplitudes to expect
 */
BOOST_AUTO_TEST_CASE(test_amplitude)
{
    // Read root files to find decay parameters
    boost::filesystem::path currentDir = boost::filesystem::path(__FILE__).parent_path();
    boost::filesystem::path cfPath("dBarCf.root");
    TFile*                  cfFile = new TFile((currentDir / cfPath).string().c_str());

    boost::filesystem::path dcsPath("dDcs.root");
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

    size_t numEvents = 1000; // Number of events in our ROOT files
    // DCS/CF magnitudes of DCS decays
    std::vector<double> dcsDcsMagnitude(numEvents);
    std::vector<double> dcsCfMagnitude(numEvents);

    // DCS/CF magnitudes of CFdecays
    std::vector<double> cfDcsMagnitude(numEvents);
    std::vector<double> cfCfMagnitude(numEvents);

    for (size_t i = 0; i < numEvents; ++i) {
        dcsDcsMagnitude[i] = std::abs(amplitude(dcs.events[i], dcsFunc));
        dcsCfMagnitude[i]  = std::abs(amplitude(dcs.events[i], cfFunc));

        cfDcsMagnitude[i] = std::abs(amplitude(cf.events[i], dcsFunc));
        cfCfMagnitude[i]  = std::abs(amplitude(cf.events[i], cfFunc));
    }

    // Average magnitudes
    double dcsDcsAvg = std::accumulate(dcsDcsMagnitude.begin(), dcsDcsMagnitude.end(), 0.0) / numEvents;
    double dcsCfAvg  = std::accumulate(dcsCfMagnitude.begin(), dcsCfMagnitude.end(), 0.0) / numEvents;

    double cfDcsAvg = std::accumulate(cfDcsMagnitude.begin(), cfDcsMagnitude.end(), 0.0) / numEvents;
    double cfCfAvg  = std::accumulate(cfCfMagnitude.begin(), cfCfMagnitude.end(), 0.0) / numEvents;

    std::cout << "DCS\nDCS component\tCF component\n";
    std::cout << dcsDcsAvg << "\t\t" << dcsCfAvg << "\n\n";

    std::cout << "CF\nDCS component\t\tCF component\n";
    std::cout << cfDcsAvg << "\t\t" << cfCfAvg << "\n";

    BOOST_CHECK(dcsDcsAvg > dcsCfAvg);
    BOOST_CHECK(cfDcsAvg < cfCfAvg);

    delete cfFile;
    delete dcsFile;
}

BOOST_AUTO_TEST_SUITE_END()
