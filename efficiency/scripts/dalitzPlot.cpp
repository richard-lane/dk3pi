/*
 * Plot a Dalitz plot from a (hard-coded) ROOT file
 *
 */
#include <iostream>

#include <TFile.h>
#include <TH2D.h>
#include <boost/filesystem.hpp>

#include "ReadRoot.h"
#include "util.h"

int main()
{
    boost::filesystem::path  currentDir = boost::filesystem::path(__FILE__).parent_path();
    boost::filesystem::path  rootFile("../../coolProjectFiles/bigRootFile.root");
    std::string              treeName          = "TupleDstToD0pi_D0ToKpipipi/DecayTree";
    std::vector<std::string> branchNames       = {"Kminus", "pi1plus", "pi2plus", "pi3minus"};
    std::vector<std::string> momentumPostfixes = {"_PX", "_PY", "_PZ", "_PE"};

    // Read the relevant data from the ROOT file
    std::unique_ptr<TFile> tFile(new TFile(rootFile.string().c_str()));
    ReadRoot               RootFile(tFile.get(), treeName, branchNames, momentumPostfixes);

    // Calculate invariant masses from these events
    size_t              numEvents = RootFile.events.size();
    std::vector<double> m12(numEvents); // m(Kpi1)
    std::vector<double> m13(numEvents); // m(Kpi2)
    for (size_t i = 0; i < numEvents; ++i) {
        m12[i] =
            invariantMass(std::vector<kinematicParams_t>{RootFile.events[i].kParams, RootFile.events[i].pi1Params});
        m13[i] =
            invariantMass(std::vector<kinematicParams_t>{RootFile.events[i].kParams, RootFile.events[i].pi2Params});
    }

    // Plot them on a stupid, annoying, impossible ROOT histogram
    // Create and fill histograms for these events (Dalitz plots)
    size_t                numBins = 100;
    double                min     = 500;
    double                max     = 1800;
    std::unique_ptr<TH2D> m12m23Hist(new TH2D("m12", "m12 vs m23", numBins, min, max, numBins, min, max));
    m12m23Hist->FillN(numEvents, m12.data(), m13.data(), nullptr);

    //  Save Dalitz plots to file
    util::saveObjectToFile(m12m23Hist.get(), "m12m23.png", "COLZ");

    return 0;
}
