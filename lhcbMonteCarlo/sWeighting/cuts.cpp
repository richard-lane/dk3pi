#include <iostream>

#include <boost/progress.hpp>

#include <TFile.h>
#include <TTree.h>

#include "cuts.h"

static std::unique_ptr<TTree> getTree(TFile& f, const std::string& treeName)
{
    std::cout << "Reading " << treeName << std::endl;
    TTree* t{nullptr};
    f.GetObject(treeName.c_str(), t);
    assert(t);

    return std::unique_ptr<TTree>(t);
}

bool passesDMassCut(const double dMass)
{
    return std::fabs(dMass - 1864.83) < 24;
}

bool passesDPTCut(const double dPT)
{
    return dPT > 4940;
}

bool passesPiSoftGhostCut(const double pGhost)
{
    return pGhost < 0.05;
}

bool passesPiSoftPTCut(const double piPT)
{
    return piPT > 360;
}

void promptCuts(const std::string& inFile,
                const std::string& treeName,
                const std::string& outFile,
                const std::string& piSoftPTBranchName,
                const std::string& piSoftGhostProbBranchName,
                const std::string& dPtBranchName,
                const std::string& dMassBranchName)
{
    // Read old tree
    TFile in(inFile.c_str(), "READ");
    auto  inTree = getTree(in, treeName);

    // Set branch addresses
    double piSoftPt{};
    inTree->SetBranchAddress(piSoftPTBranchName.c_str(), &piSoftPt);

    double piSoftGhostProb{};
    inTree->SetBranchAddress(piSoftGhostProbBranchName.c_str(), &piSoftGhostProb);

    double dPT{};
    inTree->SetBranchAddress(dPtBranchName.c_str(), &dPT);

    double dMass{};
    inTree->SetBranchAddress(dMassBranchName.c_str(), &dMass);

    // Create new tree + file
    TFile                  out(outFile.c_str(), "RECREATE");
    std::unique_ptr<TTree> outTree(inTree->CloneTree(0)); // Don't copy events

    // Lambda for our cuts
    auto passes_cuts = [&](const double piSoftPt, const double piSoftGhostProb, const double dPT, const double dMass) {
        return passesDMassCut(dMass) && passesDPTCut(dPT) && passesPiSoftGhostCut(piSoftGhostProb) &&
               passesPiSoftPTCut(piSoftPt);
    };

    // Iterate over events checking if cuts passed
    auto                    nEntries = inTree->GetEntries();
    boost::progress_display progressBar(nEntries);
    for (decltype(nEntries) i{0}; i < nEntries; ++i) {
        inTree->GetEntry(i);
        if (passes_cuts(piSoftPt, piSoftGhostProb, dPT, dMass)) {
            outTree->Fill();
        }
        ++progressBar;
    }

    // Write new tree
    out.Write();
}

void semileptonicCuts(const std::string& inFile,
                      const std::string& treeName,
                      const std::string& outFile,
                      const std::string& dPtBranchName,
                      const std::string& dMassBranchName)
{
    // Read old tree
    TFile in(inFile.c_str(), "READ");
    auto  inTree = getTree(in, treeName);

    // Set branch addresses
    double dPT{};
    inTree->SetBranchAddress(dPtBranchName.c_str(), &dPT);

    double dMass{};
    inTree->SetBranchAddress(dMassBranchName.c_str(), &dMass);

    // Create new tree + file
    TFile                  out(outFile.c_str(), "RECREATE");
    std::unique_ptr<TTree> outTree(inTree->CloneTree(0)); // Don't copy events

    // Lambda for our cuts
    auto passes_cuts = [&](const double dPT, const double dMass) { return passesDMassCut(dMass) && passesDPTCut(dPT); };

    // Iterate over events checking if cuts passed
    auto                    nEntries = inTree->GetEntries();
    boost::progress_display progressBar(nEntries);
    for (decltype(nEntries) i{0}; i < nEntries; ++i) {
        inTree->GetEntry(i);
        if (passes_cuts(dPT, dMass)) {
            outTree->Fill();
        }
        ++progressBar;
    }

    // Write new tree
    out.Write();
}
