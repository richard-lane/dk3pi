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

#include "../include/DataSetsRatio.h"
#include "../include/DecaysData.h"
#include "DataSetsRatio.cpp"
#include "DecaysData.cpp"
#include "plottingHelpers.cpp"

/*
 * Bin the CF and Mixed decays modelled in an AmpGen generated inputFile into phase bins as defined by $BIN_LIMITS
 *
 * This function is far too long but i think its probably ok for now
 *
 */
void bin_generated_decays(TFile *mixedDecays, TFile *favouredDecays)
{

    // Create D2K3piData class instances to encapsulate the data in our ROOT file
    D2K3PiData MixedData    = D2K3PiData(mixedDecays, "DalitzEventList");
    D2K3PiData FavouredData = D2K3PiData(favouredDecays, "DalitzEventList");

    // Bin decays in phase space; find a distribution of decay times in each of these bins
    MixedData.performBinning("D_decayTime");
    FavouredData.performBinning("Dbar0_decayTime");

    // Make some plots to check that the data from ROOT has been read in correctly
    // plot_things(MixedData.kVectors, MixedData.pi1Vectors, MixedData.pi2Vectors);
    // plot_things(FavouredData.kVectors, FavouredData.pi1Vectors, FavouredData.pi2Vectors);

    // Plot the hist of times in each bin
    // for (size_t bin = 0; bin < NUM_BINS; bin++) {
    //    MixedData.plotBinnedTimes(bin);
    //    FavouredData.plotBinnedTimes(bin);
    //}

    // Create a vector of DataSetsRatio objects to compare our two datasets
    std::vector<DataSetsRatio> dataSetRatios{};
    for (int bin = 0; bin < NUM_BINS; bin++) {
        dataSetRatios.push_back(DataSetsRatio(MixedData.numPointsPerTimeBin[bin],
                                              FavouredData.numPointsPerTimeBin[bin],
                                              MixedData.binCentres,
                                              MixedData.binErrors));

        // Boolean arg tells us whether to draw the graphs (useful for debugging) or to just fit the data
        std::string title = "Phase Space Bin " + std::to_string(bin);
        dataSetRatios[bin].fitToData(true, title);
    }
}

int main(int argc, char *argv[])
{
    if (argc != 3) {
        std::cerr << "Usage: <executable> <Mixed decays ROOT file> <Favoured decays ROOT file>" << std::endl;
        throw;
    }
    TFile *mixedDecays    = new TFile(argv[1]);
    TFile *favouredDecays = new TFile(argv[2]);
    bin_generated_decays(mixedDecays, favouredDecays);
}
