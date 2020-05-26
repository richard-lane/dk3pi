/*
 * bin_generated_decays.cpp
 * ROOT macro to bin AmpGen generated D -> K3pi decays into predefined bins
 *
 * NOTE: if this fails to build with error "dcs.so: No such file or directory", try restarting ROOT and building again.
 *
 * This script contains no error handling, arg verification or anything to make it work nicely
 */

#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TRandom.h"

#include "PhaseSpaceBinning.h"
#include "RatioCalculator.h"
#include "ReadAmpGen.h"
#include "RootFitter.h"
#include "util.h"

/*
 * Bin the CF and Mixed decays modelled in an AmpGen generated inputFile into phase bins as defined by $BIN_LIMITS
 * Using new APIs
 *
 * This should be split into functions
 *
 */
void bin_generated_decays(TFile *mixedDecays, TFile *favouredDecays)
{
    // Create D2K3PiData classes to represent our data and populate them.
    D2K3PiData MixedData    = D2K3PiData(mixedDecays, "DalitzEventList");
    D2K3PiData FavouredData = D2K3PiData(favouredDecays, "DalitzEventList");
    std::cout << "Reading D ROOT file" << std::endl;
    MixedData.populate("D_decayTime");
    std::cout << "Reading Dbar ROOT file" << std::endl;
    FavouredData.populate("Dbar0_decayTime");

    // Pass the times and kinematic data into the phase space binner and perform binning
    PhaseSpaceBinning MixedDataBinner = PhaseSpaceBinning(
        MixedData.decayTimes, MixedData.kVectors, MixedData.pi1Vectors, MixedData.pi2Vectors, MixedData.pi3Vectors);
    PhaseSpaceBinning FavouredDataBinner = PhaseSpaceBinning(FavouredData.decayTimes,
                                                             FavouredData.kVectors,
                                                             FavouredData.pi1Vectors,
                                                             FavouredData.pi2Vectors,
                                                             FavouredData.pi3Vectors);
    std::cout << "Binning D decays" << std::endl;
    MixedDataBinner.performBinning();

    std::cout << "Binning Dbar decays" << std::endl;
    FavouredDataBinner.performBinning();

    // Define time-space bin limits to perform ratio calculation with
    std::vector<double> timeBinLimits{};
    timeBinLimits = util::exponentialBinLimits(0.006, 2500, 15);

    std::cout << "Taking ratios" << std::endl;
    // Calculate the ratios between the datasets in each bin
    std::vector<RatioCalculator> ratios{};
    for (int bin = 0; bin < NUM_BINS; ++bin) {
        std::vector<size_t> cfCounts  = util::binVector(FavouredDataBinner._binnedTimes[bin], timeBinLimits);
        std::vector<size_t> dcsCounts = util::binVector(MixedDataBinner._binnedTimes[bin], timeBinLimits);
        ratios.push_back(RatioCalculator(cfCounts, dcsCounts, timeBinLimits));
        ratios[bin].calculateRatios();
    }

    // Fit each ratio to a plot
    std::cout << "performing fits" << std::endl;
    std::vector<RootFitter> fitters{};
    for (int bin = 0; bin < NUM_BINS; ++bin) {
        FitData_t thisBinFitData = FitData(timeBinLimits, ratios[bin].ratio, ratios[bin].error);
        fitters.push_back(RootFitter(thisBinFitData));
        fitters[bin].fit(0, 1); // Might need to change the max time here
        std::string title = "PhaseSpaceBin" + std::to_string(bin) + ".pdf";
        fitters[bin].saveFitPlot("Phase Space Bin " + std::to_string(bin), title);
    }
}

// Hide this program's main() function from ROOT's Cling interpreter
#ifndef __CINT__
int main(int argc, char *argv[])
{
    if (argc != 3) {
        std::cerr << "Usage: <executable> <Mixed decays ROOT file> <Favoured decays ROOT file>" << std::endl;
        throw;
    }
    TFile *mixedDecays    = new TFile(argv[1]);
    TFile *favouredDecays = new TFile(argv[2]);
    bin_generated_decays(mixedDecays, favouredDecays);

    return 0;
}
#endif // __CINT__
