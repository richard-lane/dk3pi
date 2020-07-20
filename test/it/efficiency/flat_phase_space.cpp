#define BOOST_TEST_DYN_LINK
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <random>

#include <TFile.h>
#include <TGenPhaseSpace.h>
#include <TH2D.h>

#include "ReadRoot.h"
#include "amplitudes.h"
#include "efficiencyUtil.h"
#include "flatPhsp.h"
#include "util.h"

BOOST_AUTO_TEST_SUITE(test_phase_space)

/*
 * Check that a randomly generated 2-body decay (flat in phsp) is produced with the correct invariant mass
 */
BOOST_AUTO_TEST_CASE(test_two_body_decay)
{
    // Decay is D -> Kpi
    TLorentzVector        dMomentum(0.0, 0.0, 0.0, D_MASS_GEV);
    std::array<double, 2> finalStateMasses = {K_MASS_GEV, PI_MASS_GEV};
    TGenPhaseSpace        PhaseSpace;
    PhaseSpace.SetDecay(dMomentum, finalStateMasses.size(), finalStateMasses.data());

    // Random number stuff
    std::random_device                     rd;
    std::mt19937                           gen(rd());
    std::uniform_real_distribution<double> uniformDistribution(0, PhaseSpace.GetWtMax());

    // Generate 5 events, check all of them have the right invariant mass
    for (size_t i = 0; i < 5; ++i) {
        std::vector<kinematicParams_t> event = randomEvent(PhaseSpace, &gen, uniformDistribution);
        BOOST_CHECK(invariantMass(event) == D_MASS_GEV);
    }
}

/*
 * Return a pointer to the amplitude model in the provided location, relative to this file
 */
amplitudeFcnPtr amplitudeModel(const std::string& path)
{
    boost::filesystem::path currentDir = boost::filesystem::path(__FILE__).parent_path();
    boost::filesystem::path lib(path);

    return readFromSharedLib((currentDir / lib).string());
}

/*
 * Read in a D->K3pi ROOT file and return a class instance containing the relevant information
 */
ReadRoot readAmpGenFile(const std::string& path)
{
    boost::filesystem::path currentDir = boost::filesystem::path(__FILE__).parent_path();
    boost::filesystem::path relPath(path);
    TFile*                  tfile = new TFile((currentDir / relPath).c_str());

    // These are the branch names that ampgen chooses for K+3pi-
    std::vector<std::string> branchNames = {"_1_K~", "_2_pi#", "_3_pi#", "_4_pi~"};

    ReadRoot data(tfile, "DalitzEventList", branchNames);

    delete tfile;
    return data;
}

/*
 * Generate event uniformly-distributed in phase space, modulate by AmpGen's amplitude model and then compare them to
 * events generated with AmpGen
 */
BOOST_AUTO_TEST_CASE(test_phase_space)
{
    // Read in AmpGen-generated DCS and CF ROOT files and store the relevant stuff in a class instance
    ReadRoot cf  = readAmpGenFile("dBarCf.root");
    ReadRoot dcs = readAmpGenFile("dDcs.root");

    // Find DCS and CF amplitude functions
    // These are just pointers to functions that work out our amplitude at a given phase space point
    amplitudeFcnPtr dcsFunc = amplitudeModel("../../../AmpGenTools/amplitude_models/dcs.so");
    amplitudeFcnPtr cfFunc  = amplitudeModel("../../../AmpGenTools/amplitude_models/cf.so");

    // Count how many event there are
    const size_t numEvents = cf.events.size();
    if (numEvents != dcs.events.size()) {
        BOOST_CHECK_MESSAGE(
            false,
            "Must have the same number of DCS and CF events for flat phase space test (just to make things nicer)");
    }

    // Generate events that should be flat in phase space
    std::random_device    rd;
    std::mt19937          gen(rd());
    std::vector<dDecay_t> flatEvents = flatDk3pi(numEvents, &gen);

    // Find the values of the DCS and CF amplitudes at these phase-space points
    std::vector<double> dcsMatrixElement(numEvents);
    std::vector<double> cfMatrixElement(numEvents);
    for (size_t i = 0; i < numEvents; ++i) {
        dcsMatrixElement[i] = std::norm(amplitude(flatEvents[i], dcsFunc));
        cfMatrixElement[i]  = std::norm(amplitude(flatEvents[i], cfFunc));
    }

    // Scale these amplitudes to have average values of 1, so that we don't affect the total number of events by
    // weighting
    double dcsAvg = std::accumulate(dcsMatrixElement.begin(), dcsMatrixElement.end(), 0.0) / dcsMatrixElement.size();
    double cfAvg  = std::accumulate(cfMatrixElement.begin(), cfMatrixElement.end(), 0.0) / cfMatrixElement.size();
    std::transform(dcsMatrixElement.begin(),
                   dcsMatrixElement.end(),
                   dcsMatrixElement.begin(),
                   std::bind(std::divides<double>(), std::placeholders::_1, dcsAvg));

    std::transform(cfMatrixElement.begin(),
                   cfMatrixElement.end(),
                   cfMatrixElement.begin(),
                   std::bind(std::divides<double>(), std::placeholders::_1, cfAvg));

    // Choose some bins in our phase space parameters to use
    const size_t nBins = 10;
    const double minX  = 0.2;
    const double minY  = 0.2;
    const double maxX  = 1.4;
    const double maxY  = 1.4;

    // 2D histograms to store our phase space projections
    TH2D* dcsHistogram =
        new TH2D("DCS Flat Events", "DCS \"Flat Events\" with weights", nBins, minX, maxX, nBins, minY, maxY);
    TH2D* cfHistogram =
        new TH2D("CF Flat Events", "CF \"Flat Events\" with weights", nBins, minX, maxX, nBins, minY, maxY);

    // Calculate the relevant invariant masses and fill in the histogram with the weights calculated from the amplitude
    // model
    for (size_t i = 0; i < numEvents; ++i) {
        // Weight the flat events with the amplitudes^2 calculated from the amplitude models
        dcsHistogram->Fill(
            invariantMass(std::vector<kinematicParams_t>{flatEvents[i].pi1Params, flatEvents[i].pi3Params}),
            invariantMass(std::vector<kinematicParams_t>{flatEvents[i].pi3Params, flatEvents[i].pi2Params}),
            dcsMatrixElement[i]);

        cfHistogram->Fill(
            invariantMass(std::vector<kinematicParams_t>{flatEvents[i].pi1Params, flatEvents[i].pi3Params}),
            invariantMass(std::vector<kinematicParams_t>{flatEvents[i].pi3Params, flatEvents[i].pi2Params}),
            cfMatrixElement[i]);
    }

    // Creat a histogram with the same bins for AmpGen events
    TH2D* ampgenDcsHist = new TH2D("DCS Ampgen Events", "DCS AmpGen", nBins, minX, maxX, nBins, minY, maxY);
    TH2D* ampgenCfHist  = new TH2D("CF Ampgen Events", "CF AmpGen", nBins, minX, maxX, nBins, minY, maxY);
    for (size_t i = 0; i < numEvents; ++i) {
        // This time we don't need to weight our events
        ampgenDcsHist->Fill(
            invariantMass(std::vector<kinematicParams_t>{dcs.events[i].pi1Params, dcs.events[i].pi3Params}),
            invariantMass(std::vector<kinematicParams_t>{dcs.events[i].pi3Params, dcs.events[i].pi2Params}));

        ampgenCfHist->Fill(
            invariantMass(std::vector<kinematicParams_t>{cf.events[i].pi1Params, cf.events[i].pi3Params}),
            invariantMass(std::vector<kinematicParams_t>{cf.events[i].pi3Params, cf.events[i].pi2Params}));
    }

    // Plot the histograms
    util::saveObjectToFile(dcsHistogram, "dcs.png");
    util::saveObjectToFile(cfHistogram, "cf.png");
    util::saveObjectToFile(ampgenDcsHist, "dcsAmpgen.png");
    util::saveObjectToFile(ampgenCfHist, "cfAmpgen.png");

    // Do a chisquared test to check if these histograms are the same
    std::cout << ampgenCfHist->Chi2Test(cfHistogram, "NORM P CHI2/NDF");

    // Enforce that our chi sq between CF and DCS is not below some special value
    double criticalChiSq = 1.05;
    BOOST_CHECK(ampgenCfHist->Chi2Test(dcsHistogram, "NORM P CHI2/NDF") > criticalChiSq);
    BOOST_CHECK(ampgenDcsHist->Chi2Test(cfHistogram, "NORM P CHI2/NDF") > criticalChiSq);

    // Enforce that our chisq/Ndf is below some critical value
    double cfChiSq  = ampgenCfHist->Chi2Test(cfHistogram, "NORM P CHI2/NDF");
    double dcsChiSq = ampgenDcsHist->Chi2Test(dcsHistogram, "NORM P CHI2/NDF");
    BOOST_CHECK_MESSAGE(cfChiSq < criticalChiSq, "CF ChiSq " + std::to_string(cfChiSq));
    BOOST_CHECK_MESSAGE(dcsChiSq < criticalChiSq, "DCS ChiSq " + std::to_string(dcsChiSq));

    delete cfHistogram;
    delete dcsHistogram;
    delete ampgenDcsHist;
    delete ampgenCfHist;
}

BOOST_AUTO_TEST_SUITE_END()
