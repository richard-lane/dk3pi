/*
 * Plot a Dalitz plot from a (hard-coded) ROOT file
 *
 */
#include <iostream>
#include <utility>

#include <TFile.h>
#include <TH2D.h>
#include <boost/filesystem.hpp>

#include "ReadRoot.h"
#include "util.h"

/*
 * Plot a Dalitz projection of a D->K3pi ROOT file
 */
void dalitzPlot(const std::string&               rootFilePath,
                const std::string&               treeName,
                const std::vector<std::string>&  branchNames,
                const std::vector<std::string>&  momentumPostfixes,
                const std::string&               plotPath,
                const std::pair<double, double>& axisLimits = std::make_pair<double, double>(500., 1800.))
{
    // Read the relevant data from the ROOT file
    std::unique_ptr<TFile> tFile(new TFile(rootFilePath.c_str()));
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
    size_t                numBins = 100;
    double                min     = axisLimits.first;
    double                max     = axisLimits.second;
    std::unique_ptr<TH2D> m12m23Hist(new TH2D(
        "dalitz", std::string(rootFilePath + " m12 vs m23;m12;m23").c_str(), numBins, min, max, numBins, min, max));
    m12m23Hist->FillN(numEvents, m12.data(), m13.data(), nullptr);

    //  Save Dalitz plots to file
    util::saveObjectToFile(m12m23Hist.get(), plotPath.c_str(), "COLZ");
}

int main()
{
    // AmpGen
    std::string              ampgenRootFile("../../AmpGen/binning/Mixed.root");
    std::string              ampgenTreeName          = "DalitzEventList";
    std::vector<std::string> ampgenBranchNames       = {"_1_K~", "_2_pi#", "_3_pi#", "_4_pi~"};
    std::vector<std::string> ampgenMomentumPostfixes = {"_Px", "_Py", "_Pz", "_E"};

    dalitzPlot(ampgenRootFile,
               ampgenTreeName,
               ampgenBranchNames,
               ampgenMomentumPostfixes,
               "ampgen.png",
               std::make_pair<double, double>(0.5, 1.8));

    // LHCb Monte Carlo
    std::string              mcRootFile("../../coolProjectFiles/bigRootFile.root");
    std::string              mcTreeName          = "TupleDstToD0pi_D0ToKpipipi/DecayTree";
    std::vector<std::string> mcBranchNames       = {"Kminus", "pi1plus", "pi2plus", "pi3minus"};
    std::vector<std::string> mcMomentumPostfixes = {"_PX", "_PY", "_PZ", "_PE"};

    dalitzPlot(mcRootFile, mcTreeName, mcBranchNames, mcMomentumPostfixes, "mc.png");

    // Real data?
    std::string              realRootFile("../../coolProjectFiles/00083875_00000123_1.charm_d02hhhh_dvntuple.root");
    std::string              realTreeName          = "Hlt2Dstp2D0Pip_D02KpPimPimPip_Tuple/DecayTree";
    std::vector<std::string> realBranchNames       = {"D0_P0", "D0_P1", "D0_P2", "D0_P3"};
    std::vector<std::string> realMomentumPostfixes = {"_PX", "_PY", "_PZ", "_PE"};

    dalitzPlot(realRootFile, realTreeName, realBranchNames, realMomentumPostfixes, "real.png");

    return 0;
}
