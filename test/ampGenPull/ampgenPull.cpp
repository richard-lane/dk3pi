#include "D2K3PiError.h"
#include "FitterUtils.h"
#include "PhysicalFitter.h"
#include "RatioCalculator.h"
#include "ReadAmpGen.h"
#include "util.h"

#include <TFile.h>

#include <boost/filesystem.hpp>

#include <fstream>
#include <memory>

/*
 * Read in two ROOT files d.root and dBar.root, take their ratios of decay times
 */
void ampgenFit(const char* dFile, const char* dBarFile, const char* outFile)
{
    // Data that will be used for fitting/calculating pulls
    double width = 2438.4;
    double x     = 0.0037;
    double y     = 0.0066;

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
    FitData_t         fitData(binLimits, ratios.ratio, ratios.error);
    IntegralOptions_t integralOptions(true, width, binLimits, 0);
    PhysicalFitter    fitter(fitData, integralOptions, false);

    std::vector<double> initialParamGuess{x, y, 0.5, 0.5, 0.5, width};
    std::vector<double> initialErrorsGuess{1, 1, 1, 1, 1, 1};
    fitter.setFitParams(initialParamGuess, initialErrorsGuess);
    fitter.fixParameters(std::vector<std::string>{"width", "x", "y"});
    fitter.fit();

    // Save the fit parameters and errors to a text file
    std::string r    = std::to_string(fitter.fitParams.fitParams[2]);
    std::string dr   = std::to_string(fitter.fitParams.fitParamErrors[2]);
    std::string reZ  = std::to_string(fitter.fitParams.fitParams[3]);
    std::string dreZ = std::to_string(fitter.fitParams.fitParamErrors[3]);
    std::string imZ  = std::to_string(fitter.fitParams.fitParams[4]);
    std::string dimZ = std::to_string(fitter.fitParams.fitParamErrors[4]);

    std::string   resultStr = r + "," + dr + "," + reZ + "," + dreZ + "," + imZ + "," + dimZ;
    std::ofstream resultsFile;
    resultsFile.open(outFile, std::ios_base::app);
    resultsFile << resultStr << std::endl;
}

int main(int argc, char* argv[])
{
    // Should probably also check that D and Dbar have been passed in the right order
    // but what are the chances of that going wrong
    if (argc != 4) {
        std::cerr << "Usage: ./ampgenpull <d root file> <d bar root file> <outfile>" << std::endl;
        throw D2K3PiException();
    }

    // Check that the file we'll write to exists
    if (!boost::filesystem::exists(argv[3])) {
        std::cerr << argv[3] << " does not exist; will be unable to append fit results to it" << std::endl;
        throw D2K3PiException();
    }

    ampgenFit(argv[1], argv[2], argv[3]);
    return 0;
}
