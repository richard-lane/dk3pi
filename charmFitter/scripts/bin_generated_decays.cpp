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

#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TRandom.h"

#include "PhaseSpaceBinning.h"
#include "RatioCalculator.h"
#include "ReadAmpGen.h"
#include "RootFitter.h"
#include "fitterUtil.h"
#include "util.h"

/*
 * Bin the CF and Mixed decays modelled in an AmpGen generated inputFile into phase bins as defined by BIN_LIMITS
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

    // Put our data into a stupid struct
    const size_t    numBins{5};
    PhspBinningData MixedPhspData(
        MixedData.decayTimes, MixedData.kVectors, MixedData.pi1Vectors, MixedData.pi2Vectors, MixedData.pi3Vectors);
    PhspBinningData FavouredPhspData(FavouredData.decayTimes,
                                     FavouredData.kVectors,
                                     FavouredData.pi1Vectors,
                                     FavouredData.pi2Vectors,
                                     FavouredData.pi3Vectors);

    // Pass the times and kinematic data into the phase space binner and perform binning
    // These are sensitive to where the script is run from... which i could fix but i won't
    std::string       dcsModel{"../AmpGenTools/amplitude_models/dcs.so"};
    std::string       cfModel{"../AmpGenTools/amplitude_models/cf.so"};
    PhaseSpaceBinning MixedDataBinner    = PhaseSpaceBinning(MixedPhspData, dcsModel, cfModel);
    PhaseSpaceBinning FavouredDataBinner = PhaseSpaceBinning(FavouredPhspData, dcsModel, cfModel);

    std::cout << "Binning D decays" << std::endl;
    std::vector<PhspBinningData> BinnedMixedData = MixedDataBinner.performBinning();

    std::cout << "Binning Dbar decays" << std::endl;
    std::vector<PhspBinningData> BinnedFavouredData = FavouredDataBinner.performBinning();

    // Define time-space bin limits to perform ratio calculation with
    std::vector<double> timeBinLimits = FitterUtil::exponentialBinLimits(0.006, 2500, 5);

    std::cout << "Taking ratios" << std::endl;
    // Calculate the ratios between the datasets in each bin
    std::vector<std::pair<std::vector<double>, std::vector<double>>> ratiosAndErrors{};
    for (size_t bin = 0; bin < numBins; ++bin) {
        std::vector<size_t> cfCounts  = util::binVector(BinnedFavouredData[bin].decayTimes, timeBinLimits);
        std::vector<size_t> dcsCounts = util::binVector(BinnedMixedData[bin].decayTimes, timeBinLimits);
        ratiosAndErrors.push_back(RatioCalculator::ratioAndError(cfCounts, dcsCounts));
    }

    // Fit each ratio to a plot
    std::cout << "performing fits" << std::endl;
    std::vector<RootFitter> fitters{};
    for (size_t bin = 0; bin < numBins; ++bin) {
        FitData_t thisBinFitData = FitData(timeBinLimits, ratiosAndErrors[bin].first, ratiosAndErrors[bin].second);
        fitters.push_back(RootFitter(thisBinFitData));
        fitters[bin].fit(0, 1); // Might need to change the max time here
        std::string title = "PhaseSpaceBin" + std::to_string(bin) + ".png";
        fitters[bin].saveFitPlot("Phase Space Bin " + std::to_string(bin), title);
    }

    // Find the invariant mass of each mixed event and plot them on a scatter plot
    std::vector<std::unique_ptr<TH2D>> hists(numBins);
    for (size_t bin = 0; bin < numBins; ++bin) {
        hists[bin] =
            std::make_unique<TH2D>(std::to_string(bin).c_str(), std::to_string(bin).c_str(), 100, 0, 3, 100, 0, 3);
        for (size_t i = 0; i < BinnedMixedData[bin].decayTimes.size(); ++i) {
            // Find the invariant masses of each event
            double x = (BinnedMixedData[bin].kVectors[i] + BinnedMixedData[bin].pi1Vectors[i]).M2();
            double y = (BinnedMixedData[bin].pi1Vectors[i] + BinnedMixedData[bin].pi2Vectors[i] +
                        BinnedMixedData[bin].pi3Vectors[i])
                           .M2();
            // Add them to the bin'th histogram
            hists[bin]->Fill(x, y);
        }
    }

    std::unique_ptr<TCanvas> c = std::make_unique<TCanvas>();
    c->cd();

    std::vector<EColor> colours{kRed, kGreen, kBlue, kYellow, kBlack};
    for (size_t i = 0; i < numBins; ++i) {
        hists[i]->SetMarkerColor(colours[i]);
        hists[i]->Draw("SAME");
    }
    c->SaveAs("plot.png");
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
