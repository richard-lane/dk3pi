#include "D2K3PiError.h"
#include "RatioCalculator.h"
#include "ReadAmpGen.h"
#include "util.h"

#include <TFile.h>

#include <memory>

/*
 * Read in two ROOT files d.root and dBar.root, take their ratios of decay times
 */
void ampgenFit(const char* dFile, const char* dBarFile)
{
    // Data that will be used for fitting/calculating pulls
    double width = 2438.4;

    // Read in ROOT files
    std::unique_ptr<TFile> mixedDecays    = std::make_unique<TFile>(dFile);
    std::unique_ptr<TFile> favouredDecays = std::make_unique<TFile>(dBarFile);

    D2K3PiData MixedData    = D2K3PiData(mixedDecays.get(), "DalitzEventList");
    D2K3PiData FavouredData = D2K3PiData(favouredDecays.get(), "DalitzEventList");
    MixedData.populate("D_decayTime");
    FavouredData.populate("Dbar0_decayTime");

    // Time binning
    std::vector<double> binLimits = util::exponentialBinLimits(0.006, width, 25);
    std::vector<size_t> cfCounts  = util::binVector(FavouredData.decayTimes, binLimits);
    std::vector<size_t> dcsCounts = util::binVector(MixedData.decayTimes, binLimits);

    // Take ratios
    RatioCalculator ratios(cfCounts, dcsCounts, binLimits);
    ratios.calculateRatios();

    // Perform a fit

    // Save the fit parameters and errors to a text file
}

int main(int argc, char* argv[])
{
    // Should probably also check that D and Dbar have been passed in the right order
    // but what are the chances of that
    if (argc != 3) {
        std::cerr << "Usage: ./ampgenpull <d root file> <d bar root file>" << std::endl;
        throw D2K3PiException();
    }

    ampgenFit(argv[1], argv[2]);
    return 0;
}
