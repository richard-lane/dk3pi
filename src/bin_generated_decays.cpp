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

#include "../lib/Fitter.h"
#include "../lib/PhaseSpaceBinning.h"
#include "../lib/RatioCalculator.h"
#include "../lib/ReadAmpGen.h"

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
    MixedData.populate("D_decayTime");
    FavouredData.populate("Dbar0_decayTime");

    // Pass the times and kinematic data into the phase space binner and perform binning
    PhaseSpaceBinning MixedDataBinner = PhaseSpaceBinning(
        MixedData.decayTimes, MixedData.kVectors, MixedData.pi1Vectors, MixedData.pi2Vectors, MixedData.pi3Vectors);
    PhaseSpaceBinning FavouredDataBinner = PhaseSpaceBinning(FavouredData.decayTimes,
                                                             FavouredData.kVectors,
                                                             FavouredData.pi1Vectors,
                                                             FavouredData.pi2Vectors,
                                                             FavouredData.pi3Vectors);
    MixedDataBinner.performBinning();
    FavouredDataBinner.performBinning();

    // Define time-space bin limits to perform ratio calculation with
    std::vector<double> timeBinLimits{};
    for (size_t i = 0; i < 200; ++i) {
        timeBinLimits.push_back(i * i * i * 0.006 / (200 * 200 * 200));
    }

    // Calculate the ratios between the datasets in each bin
    std::vector<RatioCalculator> ratios{};
    for (int bin = 0; bin < NUM_BINS; ++bin) {
        ratios.push_back(
            RatioCalculator(FavouredDataBinner._binnedTimes[bin], MixedDataBinner._binnedTimes[bin], timeBinLimits));
        ratios[bin].calculateRatios();
    }

    // Fit each ratio to a plot
    std::vector<Fitter> fitters{};
    for (int bin = 0; bin < NUM_BINS; ++bin) {
        FitData_t thisBinFitData =
            FitData(ratios[bin].binCentres, ratios[bin].binWidths, ratios[bin].ratio, ratios[bin].error);
        fitters.push_back(Fitter(thisBinFitData));
        fitters[bin].fitUsingRootBuiltinPol2();
        std::string title = "PhaseSpaceBin" + std::to_string(bin) + ".pdf";
        fitters[bin].saveFitPlot(title, title);
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
