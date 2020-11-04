#include <experimental/filesystem>
#include <memory>
#include <regex>

#include "D2K3PiError.h"
#include "sWeighting.h"
#include "util.h"

#include <RooDstD0BG.h>
#include <RooGaussian.h>
#include <RooJohnson.h>
#include <RooPolynomial.h>
#include <RooRealVar.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TTree.h>

static void sPlotHist(const std::string&               rootFile,
                      const std::string&               plotPath,
                      const std::string&               treeName,
                      const std::string&               observableName,
                      const std::pair<double, double>& axisLimits)
{
    // Read tree
    TFile* newFile = new TFile(rootFile.c_str(), "READ");
    TTree* tree{nullptr};
    newFile->GetObject(treeName.c_str(), tree);
    assert(tree);

    // Set up buffers for reading data
    double signalWt{0};
    double backgroundWt{0};
    double observable{0};
    tree->SetBranchAddress(observableName.c_str(), &observable);
    tree->SetBranchAddress("numSignalEvents_sw", &signalWt);
    tree->SetBranchAddress("numBackgroundEvents_sw", &backgroundWt);

    // Create + fill hists
    TH1D* signal =
        new TH1D("signal", "Reweighted Events;Delta_M /MeV; Weighted Count", 100, axisLimits.first, axisLimits.second);
    TH1D* background =
        new TH1D("bkg", "Reweighted Events;Delta_M /MeV; Weighted Count", 100, axisLimits.first, axisLimits.second);
    auto numEvents{tree->GetEntries()};
    for (decltype(numEvents) i{0}; i < numEvents; ++i) {
        tree->GetEntry(i);
        signal->Fill(observable, signalWt);
        background->Fill(observable, backgroundWt);
    }

    signal->SetLineColor(kBlue);
    background->SetLineColor(kRed);
    signal->SetStats(false);

    // Legend
    util::LegendParams_t legend{0.7, 0.9, 0.7, 0.9};
    util::saveObjectsToFile<TH1D>(
        {signal, background}, {"HIST C SAME", "HIST C SAME"}, {"Signal", "Background"}, plotPath, legend);

    delete signal;
    delete background;
    delete newFile; // Can't delete this too early otherwise will segfault when reading the data
}

/*
 * Extract name of a ROOT file from its path
 *
 * Regex matchs a filename ending in .root, which may look like e.g, "file.root" "./file.root",
 * "directory/file.root", etc.
 * Should definitely be unit tested, but it isn't because i don't want to
 */
static std::string rootFileName(const std::string& path)
{
    std::smatch matches;
    std::regex_search(path, matches, std::regex(R"(.*?\/?([\w\-. ]+\.root))"));
    if (matches.empty()) {
        throw D2K3PiException();
    }

    return matches[1].str();
}

/*
 * return a vector of (branch name, units, limits) structs for branches common to both prompt and SL root files
 */
static std::vector<sWeighting::RootBranch> commonBranchParams(void)
{
    double max{4000000000};
    return {sWeighting::RootBranch{"D_IP_OWNPV", "", 0, max},   sWeighting::RootBranch{"D_IPCHI2_OWNPV", "", -max, max},
            sWeighting::RootBranch{"D_P", "MeV", 0, max},       sWeighting::RootBranch{"D_PT", "MeV", 0, max},
            sWeighting::RootBranch{"D_PE", "MeV", 0, max},      sWeighting::RootBranch{"D_PX", "MeV", -max, max},
            sWeighting::RootBranch{"D_PY", "MeV", -max, max},   sWeighting::RootBranch{"D_PZ", "MeV", -max, max},
            sWeighting::RootBranch{"D_M", "MeV", -max, max},    sWeighting::RootBranch{"D_MM", "MeV", 0, max},
            sWeighting::RootBranch{"D_MMERR", "MeV", 0, max},   sWeighting::RootBranch{"D_ID", "", -max, max},
            sWeighting::RootBranch{"D_TAU", "", -max, max},     sWeighting::RootBranch{"D_TAUERR", "", -max, max},
            sWeighting::RootBranch{"D_TAUCHI2", "", -max, max},

            sWeighting::RootBranch{"K_IP_OWNPV", "", 0, max},   sWeighting::RootBranch{"K_IPCHI2_OWNPV", "", 0, max},
            sWeighting::RootBranch{"K_P", "MeV", 0, max},       sWeighting::RootBranch{"K_PT", "MeV", 0, max},
            sWeighting::RootBranch{"K_PE", "MeV", 0, max},      sWeighting::RootBranch{"K_PX", "MeV", -max, max},
            sWeighting::RootBranch{"K_PY", "MeV", -max, max},   sWeighting::RootBranch{"K_PZ", "MeV", -max, max},
            sWeighting::RootBranch{"K_M", "MeV", -max, max},    sWeighting::RootBranch{"K_ID", "", -max, max},

            sWeighting::RootBranch{"pi1_IP_OWNPV", "", 0, max}, sWeighting::RootBranch{"pi1_IPCHI2_OWNPV", "", 0, max},
            sWeighting::RootBranch{"pi1_P", "MeV", 0, max},     sWeighting::RootBranch{"pi1_PT", "MeV", 0, max},
            sWeighting::RootBranch{"pi1_PE", "MeV", 0, max},    sWeighting::RootBranch{"pi1_PX", "MeV", -max, max},
            sWeighting::RootBranch{"pi1_PY", "MeV", -max, max}, sWeighting::RootBranch{"pi1_PZ", "MeV", -max, max},
            sWeighting::RootBranch{"pi1_M", "MeV", 0, max},     sWeighting::RootBranch{"pi1_ID", "", -max, max},

            sWeighting::RootBranch{"pi2_IP_OWNPV", "", 0, max}, sWeighting::RootBranch{"pi2_IPCHI2_OWNPV", "", 0, max},
            sWeighting::RootBranch{"pi2_P", "MeV", 0, max},     sWeighting::RootBranch{"pi2_PT", "MeV", 0, max},
            sWeighting::RootBranch{"pi2_PE", "MeV", 0, max},    sWeighting::RootBranch{"pi2_PX", "MeV", -max, max},
            sWeighting::RootBranch{"pi2_PY", "MeV", -max, max}, sWeighting::RootBranch{"pi2_PZ", "MeV", -max, max},
            sWeighting::RootBranch{"pi2_M", "MeV", 0, max},     sWeighting::RootBranch{"pi2_ID", "", -max, max},

            sWeighting::RootBranch{"pi3_IP_OWNPV", "", 0, max}, sWeighting::RootBranch{"pi3_IPCHI2_OWNPV", "", 0, max},
            sWeighting::RootBranch{"pi3_P", "MeV", 0, max},     sWeighting::RootBranch{"pi3_PT", "MeV", 0, max},
            sWeighting::RootBranch{"pi3_PE", "MeV", 0, max},    sWeighting::RootBranch{"pi3_PX", "MeV", -max, max},
            sWeighting::RootBranch{"pi3_PY", "MeV", -max, max}, sWeighting::RootBranch{"pi3_PZ", "MeV", -max, max},
            sWeighting::RootBranch{"pi3_M", "MeV", 0, max},     sWeighting::RootBranch{"pi3_ID", "", -max, max}};
}

/*
 * Return a vector of tuples (branch name, allowed range, unit) for prompt decays
 */
static std::vector<sWeighting::RootBranch> promptBranches(void)
{
    double                              max{4000000000};
    auto                                commonBranches{commonBranchParams()};
    std::vector<sWeighting::RootBranch> extraBranches{
        sWeighting::RootBranch{"Dstar_IPCHI2_OWNPV", "", -max, max},
        sWeighting::RootBranch{"Dstar_P", "MeV", 0, max},
        sWeighting::RootBranch{"Dstar_PT", "MeV", 0, max},
        sWeighting::RootBranch{"Dstar_PE", "MeV", 0, max},
        sWeighting::RootBranch{"Dstar_PX", "MeV", -max, max},
        sWeighting::RootBranch{"Dstar_PY", "MeV", -max, max},
        sWeighting::RootBranch{"Dstar_PZ", "MeV", 0, max},
        sWeighting::RootBranch{"Dstar_M", "MeV", 0, max},
        sWeighting::RootBranch{"Dstar_MM", "MeV", 0, max},
        sWeighting::RootBranch{"Dstar_ID", "", -max, max},

        sWeighting::RootBranch{"pisoft_IP_OWNPV", "", 0, max},
        sWeighting::RootBranch{"pisoft_IPCHI2_OWNPV", "", 0, max},
        sWeighting::RootBranch{"pisoft_P", "MeV", 0, max},
        sWeighting::RootBranch{"pisoft_PT", "MeV", 0, max},
        sWeighting::RootBranch{"pisoft_PE", "MeV", 0, max},
        sWeighting::RootBranch{"pisoft_PX", "MeV", -max, max},
        sWeighting::RootBranch{"pisoft_PY", "MeV", -max, max},
        sWeighting::RootBranch{"pisoft_PZ", "MeV", -max, max},
        sWeighting::RootBranch{"pisoft_M", "MeV", 0, max},
        sWeighting::RootBranch{"pisoft_ID", "", -max, max},

        sWeighting::RootBranch{"DELTA_M", "MeV", 139.3, 168.132}}; // hard coded ew

    commonBranches.insert(commonBranches.end(), extraBranches.begin(), extraBranches.end());
    return commonBranches;
}

/*
 * If no branch containing the D mass difference exists, we may have to create one
 *
 * This does that
 *
 * Copies the old root file to copy_<filename>, in case of things going bad
 */
static void addDeltaMBranch(const std::string& path,
                            const std::string& treeName,
                            const std::string& dMassBranchName,
                            const std::string& dStarMassBranchName,
                            const std::string& deltaMBranchName)
{
    std::cout << "Creating branch " << deltaMBranchName << std::endl;
    // Find the filename from the path so we know what to copy to
    std::string newPath{"copy_" + rootFileName(path)};

    // Copy old file before we do anything else, in case we break it
    std::cout << "Copying " << path << " to " << newPath << std::endl;
    std::experimental::filesystem::copy(path, newPath);

    // Read our tree
    TFile  f(path.c_str(), "update");
    TTree* tree = f.Get<TTree>(treeName.c_str());
    assert(tree);

    // Create the new branch
    double deltaM{0};
    double dMass{0};
    double dStarMass{0};
    tree->SetBranchAddress(dMassBranchName.c_str(), &dMass);
    tree->SetBranchAddress(dStarMassBranchName.c_str(), &dStarMass);
    TBranch* deltaMBranch{tree->Branch(deltaMBranchName.c_str(), &deltaM, (deltaMBranchName + "/D").c_str())};

    // Populate the branch
    auto numEvents = tree->GetEntries();
    for (decltype(numEvents) i = 0; i < numEvents; ++i) {
        tree->GetEntry(i);
        deltaM = dStarMass - dMass;
        deltaMBranch->Fill();
    }

    tree->Write("", TObject::kOverwrite);
}

/*
 * Fit to prompt D->K3Pi
 *
 * Spit a new root with with weights out at outFile; draw a plot at outImg
 */
[[maybe_unused]] static void promptFit(const std::string& inFile,
                                       const std::string& outFile,
                                       const std::string& outImg,
                                       const std::string& outMassFitPlot,
                                       const std::string& graph)
{
    const std::string treeName{"dk3pi/DecayTree"};

    // The parameters in our models that will be fixed after fitting
    const std::vector<std::string> paramsToFix{
        "a", "b", "c", "Gamma", "Delta", "Delta M0", "D Mass Width", "Mean D Mass"};

    // Create a signal model
    // These RooRealVars need to have the same scope as the models that use them, otherwise the models will segfault
    std::string  observableName{"DELTA_M"};
    const double loDeltaM{139.3};
    const double hiDeltaM{168.0};
    RooRealVar   deltaM(observableName.c_str(), observableName.c_str(), loDeltaM, hiDeltaM, "MeV/c^2");
    RooRealVar   meanDMassVar("Mean D Mass", "Mean D Mass", 146.0, 144.0, 147.0);
    RooRealVar   dMassWidthVar("D Mass Width", "D Mass Width", 1.0, 0.0001, 5.0);
    RooRealVar   gammaVar("Gamma", "Gamma", 0.01, -5.0, 5.0);
    RooRealVar   deltaVar("Delta", "Delta", 1.0, 0.0001, 10.0);
    RooJohnson   signalModel("Signal Johnson",
                           "Signal PDF: Johnson Function",
                           deltaM,
                           meanDMassVar,
                           dMassWidthVar,
                           gammaVar,
                           deltaVar,
                           135.0); // The value below which our signal PDF is set to 0

    // Create a background model
    RooRealVar deltaM0Var("Delta M0", "Delta M0", deltaM.getMin());
    RooRealVar aVar("a", "a", 0.01, -10.0, 10.0);
    RooRealVar bVar("b", "b", 0.01, -10.0, 10.0);
    RooRealVar cVar("c", "c", 2.0, 0.0001, 100.0);
    RooDstD0BG backgroundModel =
        RooDstD0BG("Prompt Background", "Background PDF: Prompt D", deltaM, deltaM0Var, aVar, bVar, cVar);

    // Add a delta_m branch to our root file
    addDeltaMBranch(inFile, treeName, "D_M", "Dstar_M", "DELTA_M");

    // sWeight the data
    // The new DELTA_M branch gets written to a tree at /DecayTree/, not /dk3pi/DecayTree/ so i guess let's use that
    std::unique_ptr<TTree> weightedTree = sWeighting::createSWeightedTree(inFile,
                                                                          "DecayTree",
                                                                          promptBranches(),
                                                                          signalModel,
                                                                          backgroundModel,
                                                                          observableName,
                                                                          paramsToFix,
                                                                          outMassFitPlot.c_str(),
                                                                          graph.c_str());

    // Write the weighted tree to a new file
    weightedTree->SaveAs(outFile.c_str());

    // Read in the new ROOT file and plot a histogram of reweighted mass difference
    sPlotHist(outFile, outImg.c_str(), "RooTreeDataStore_data_DecayTree", "DELTA_M", {loDeltaM, hiDeltaM});
}

/*
 * Fit to semileptonic D -> K3pi
 *
 * Spit a new root with with weights out at outFile; draw a plot at outImg
 */
[[maybe_unused]] static void semiLeptonicFit(const std::string& inFile,
                                             const std::string& outFile,
                                             const std::string& outImg,
                                             const std::string& outMassFitPlot,
                                             const std::string& graph)
{
    const std::string treeName{"dk3pi/DecayTree"};

    // The parameters in our models that will be fixed after fitting
    const std::vector<std::string> paramsToFix{"a", "b", "D_M", "D Mass Width", "Mean D Mass"};

    // Create a signal model
    // These RooRealVars need to have the same scope as the models that use them, otherwise the models will segfault
    std::string  observableName{"D_M"};
    const double lowDMass{1825.};
    const double hiDMass{1905.};
    RooRealVar   dMassVar(observableName.c_str(), observableName.c_str(), lowDMass, hiDMass, "MeV/c^2");
    RooRealVar   meanDMassVar("Mean D Mass", "Mean D Mass", 1865.0, 1850.0, 1880.0);
    RooRealVar   dMassWidthVar("D Mass Width", "D Mass Width", 10.0, 0.0001, 50.0);
    RooGaussian  signalModel("Signal Gaussian", "Signal PDF: Gaussian Function", dMassVar, meanDMassVar, dMassWidthVar);

    // Create a background model
    RooRealVar    aVar("a", "a", 0.005, 0.0, 0.01);
    RooRealVar    bVar("b", "b", 0.0);
    RooPolynomial backgroundModel("SL Background", "Background PDF: Polynomial", dMassVar, RooArgList(aVar, bVar));

    // Create a new ROOT file with the tree not in a directory
    // Hack for now until i get the sWeighting to properly deal with directories; TODO
    std::string newFileName{"new_" + inFile};
    TFile       oldFile(inFile.c_str());
    TTree*      oldTree{nullptr};
    oldFile.GetObject(treeName.c_str(), oldTree);
    assert(oldTree);
    TFile*                newFile = new TFile(newFileName.c_str(), "RECREATE");
    [[maybe_unused]] auto newTree{oldTree->CloneTree()}; // not actually unused. i think
    newFile->Write();
    delete oldTree;
    delete newFile;

    // sWeight the data
    std::unique_ptr<TTree> weightedTree = sWeighting::createSWeightedTree(newFileName,
                                                                          "DecayTree",
                                                                          commonBranchParams(),
                                                                          signalModel,
                                                                          backgroundModel,
                                                                          observableName,
                                                                          paramsToFix,
                                                                          outMassFitPlot.c_str(),
                                                                          graph.c_str());

    // Write the weighted tree to a new file
    weightedTree->SaveAs(outFile.c_str());

    // Read in the new ROOT file and plot a hist of reweighted D mass
    sPlotHist(outFile, outImg.c_str(), "RooTreeDataStore_data_DecayTree", "D_M", {lowDMass, hiDMass});
}

int main()
{
    // const std::string rsPath{"./test_2011_RS_prompt.root"};
    // promptFit(rsPath, "newRS.root", "rs.png", "rsMassFit.png", "rsGraph.dot");

    // const std::string wsPath{"./test_2011_WS_prompt.root"};
    // promptFit(wsPath, "newWS.root", "ws.png", "wsMassFit.png", "wsGraph.dot");

    const std::string slPath{"all.root"};
    semiLeptonicFit(slPath, "newSL.root", "sl.png", "slMassFit.png", "sl.dot");

    return 0;
}
