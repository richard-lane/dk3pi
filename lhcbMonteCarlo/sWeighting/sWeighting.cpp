/*
 * Based on Jenny Smallwood's (Oxford) code that she sent as an example of D->K3Pi sWeighting
 *
 * It's a bit hacky and hard to read at the moment but i'll fix it up maybe
 * I can probably get rid of the roo workspace thing
 */
#include "sWeighting.h"

#include <experimental/filesystem>
#include <memory>
#include <regex>
#include <stdexcept>
#include <tuple>

#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooDstD0BG.h>
#include <RooJohnson.h>
#include <RooRealVar.h>
#include <RooStats/SPlot.h>
#include <RooWorkspace.h>
#include <TCanvas.h>
#include <TFile.h>

namespace sWeighting
{

struct BranchMismatch : public std::exception {
    const char* what() const throw() { return "Number of ROOT branches, allowed ranges and units did not match."; };
};

struct BadFileName : public std::runtime_error {
    BadFileName(const std::string& file)
        : std::runtime_error(std::string("Could not find a ROOT file at path " + file).c_str())
    {
    }
};

RooJohnson johnsonSignalModel(RooAbsReal& dMassRange,
                              RooAbsReal& meanDMass,
                              RooAbsReal& dMassWidth,
                              RooAbsReal& gamma,
                              RooAbsReal& delta,
                              double      massThreshhold)
{
    return RooJohnson("Signal Johnson",
                      "Signal PDF: Johnson Function",
                      dMassRange,
                      meanDMass,
                      dMassWidth,
                      gamma,
                      delta,
                      massThreshhold);
};

RooDstD0BG promptBackgroundModel(RooAbsReal& deltaM, RooAbsReal& deltaM0, RooAbsReal& a, RooAbsReal& b, RooAbsReal& c)
{
    // Yes, the function signature has the shape parameters in the order c, a, b.
    return RooDstD0BG("Prompt background", "Background PDF: Prompt D", deltaM, deltaM0, c, a, b);
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
        throw BadFileName(path);
    }

    return matches[1].str();
}

void addDeltaMBranch(const std::string& path,
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
 * Read data from a ROOT file into a RooWorkspace
 *
 * Need to pass in vectors of branch names, ranges (i.e. min/max allowed values) and units (e.g. "MeV^2")
 */
static void readData(RooWorkspace&                                 workspace,
                     const std::string&                            rootFilePath,
                     const std::string&                            treeName,
                     const std::vector<std::string>&               branchNames,
                     const std::vector<std::pair<double, double>>& branchRanges,
                     const std::vector<std::string>&               units)
{
    std::cout << "Reading data from" << rootFilePath << std::endl;
    // Check that we have the right number of branch names and ranges
    const size_t numBranches{branchNames.size()};
    if (numBranches != branchRanges.size() || numBranches != units.size()) {
        throw BranchMismatch();
    }

    // Create a dataset holding all the branches we want
    std::unique_ptr<RooArgSet> branches = std::make_unique<RooArgSet>();
    for (size_t i = 0; i < numBranches; ++i) {
        RooRealVar branch(branchNames[i].c_str(),
                          branchNames[i].c_str(),
                          branchRanges[i].first,
                          branchRanges[i].second,
                          units[i].c_str());
        assert(branches->addClone(branch));
    }

    // Create a RooDataSet, telling it compress data into a Tree
    RooAbsData::setDefaultStorageType(RooAbsData::Tree);
    std::unique_ptr<RooDataSet> Data = std::make_unique<RooDataSet>(
        "data", "data", *branches, RooFit::ImportFromFile(rootFilePath.c_str(), treeName.c_str()));

    // Add data to the workspace
    workspace.import(*Data, RooFit::Rename("workspace data"));
}

/*
 * Perform sWeighting
 *
 * The sPlot method requires us to fix all parameters of the model that are not yields once we have performed the fit;
 * these should be passed by name via fixedParameters
 *
 * NB: doesn't check that sufficient/sensible parameters are fixed, since that's annoyingly hard
 *
 * if a plot title C-string is provided then a plot of the mass fit will be created
 *
 */
static TTree* sWeightData(RooWorkspace&                   workspace,
                          RooAbsPdf&                      signalModel,
                          RooAbsPdf&                      backgroundModel,
                          const std::vector<std::string>& fixedParameters,
                          const char*                     plotTitle = nullptr)
{
    std::cout << "Performing sWeighting" << std::endl;
    // Create variables tracking the number of signal and number of background events
    // The numbers here are just numbers that I have chosen
    // They should probably be different and maybe be changed
    RooRealVar numSignalEvents("numSignalEvents", "Signal Events", 20000, 0, 1000000);
    RooRealVar numBackgroundEvents("numBackgroundEvents", "Background Events", 20000, 0, 1000000);

    // Create a combined model
    RooAddPdf combinedModel("combined model",
                            "Combined Signal + Background Models",
                            RooArgList(signalModel, backgroundModel),
                            RooArgList(numSignalEvents, numBackgroundEvents));

    // Add the model to the workspace
    workspace.import(combinedModel);

    // Extract the data from the workspace
    RooAbsData* workspaceData = workspace.data("workspace data");
    assert(workspaceData);
    RooDataSet data = *static_cast<RooDataSet*>(workspaceData);

    // Create a graphviz .dot file to visualise the fitter, which might be useful
    combinedModel.graphVizTree("full.dot");

    // Fit the model to data, using an extended likelihood fit
    combinedModel.fitTo(data, RooFit::Extended());

    if (plotTitle) {
        RooRealVar* deltaM{workspace.var("DELTA_M")};
        TCanvas*    c     = new TCanvas("data", "data");
        RooPlot*    frame = new RooPlot("Mass Fit", "Mass Fit", *deltaM, 130, 170, 100);
        data.plotOn(frame);
        combinedModel.plotOn(frame);
        combinedModel.plotOn(frame, RooFit::Components(signalModel), RooFit::LineStyle(kDashed));
        combinedModel.plotOn(frame, RooFit::Components(backgroundModel), RooFit::LineStyle(kDotted));
        frame->Draw();
        c->SaveAs(plotTitle);
        delete frame;
        delete c;
    }

    // Fix the parameters that we were meant to fix
    for (const auto& paramName : fixedParameters) {
        workspace.var(paramName.c_str())->setConstant();
    }

    // Perform sPlot fit to find the number of signal and background events
    RooStats::SPlot("sData", "An sPlot", data, &combinedModel, RooArgList(numSignalEvents, numBackgroundEvents));

    // TODO find out if this is the right dataset
    return data.GetClonedTree();
}

/*
 * Return a vector of tuples (branch name, allowed range, unit)
 */
static std::vector<std::tuple<std::string, std::pair<double, double>, std::string>> desiredBranches(void)
{
    std::pair<double, double> defaultUnsignedLimits(0., 4000000000.);
    std::pair<double, double> defaultSignedLimits(-4000000000., 4000000000.);
    return {std::tuple("Dstar_IPCHI2_OWNPV", defaultSignedLimits, ""),
            std::tuple("Dstar_P", defaultUnsignedLimits, "MeV"),
            std::tuple("Dstar_PT", defaultUnsignedLimits, "MeV"),
            std::tuple("Dstar_PE", defaultUnsignedLimits, "MeV"),
            std::tuple("Dstar_PX", defaultSignedLimits, "MeV"),
            std::tuple("Dstar_PY", defaultSignedLimits, "MeV"),
            std::tuple("Dstar_PZ", defaultUnsignedLimits, "MeV"),
            std::tuple("Dstar_M", defaultUnsignedLimits, "MeV"),
            std::tuple("Dstar_MM", defaultUnsignedLimits, "MeV"),
            std::tuple("Dstar_ID", defaultSignedLimits, ""),

            std::tuple("D_IP_OWNPV", defaultUnsignedLimits, ""),
            std::tuple("D_IPCHI2_OWNPV", defaultSignedLimits, ""),
            std::tuple("D_P", defaultUnsignedLimits, "MeV"),
            std::tuple("D_PT", defaultUnsignedLimits, "MeV"),
            std::tuple("D_PE", defaultUnsignedLimits, "MeV"),
            std::tuple("D_PX", defaultSignedLimits, "MeV"),
            std::tuple("D_PY", defaultSignedLimits, "MeV"),
            std::tuple("D_PZ", defaultSignedLimits, "MeV"),
            std::tuple("D_M", defaultSignedLimits, "MeV"),
            std::tuple("D_MM", defaultUnsignedLimits, "MeV"),
            std::tuple("D_MMERR", defaultUnsignedLimits, "MeV"),
            std::tuple("D_ID", defaultSignedLimits, ""),
            std::tuple("D_TAU", defaultSignedLimits, ""),
            std::tuple("D_TAUERR", defaultSignedLimits, ""),
            std::tuple("D_TAUCHI2", defaultSignedLimits, ""),

            std::tuple("K_IP_OWNPV", defaultUnsignedLimits, ""),
            std::tuple("K_IPCHI2_OWNPV", defaultUnsignedLimits, ""),
            std::tuple("K_P", defaultUnsignedLimits, "MeV"),
            std::tuple("K_PT", defaultUnsignedLimits, "MeV"),
            std::tuple("K_PE", defaultUnsignedLimits, "MeV"),
            std::tuple("K_PX", defaultSignedLimits, "MeV"),
            std::tuple("K_PY", defaultSignedLimits, "MeV"),
            std::tuple("K_PZ", defaultSignedLimits, "MeV"),
            std::tuple("K_M", defaultSignedLimits, "MeV"),
            std::tuple("K_ID", defaultSignedLimits, ""),

            std::tuple("pi1_IP_OWNPV", defaultUnsignedLimits, ""),
            std::tuple("pi1_IPCHI2_OWNPV", defaultUnsignedLimits, ""),
            std::tuple("pi1_P", defaultUnsignedLimits, "MeV"),
            std::tuple("pi1_PT", defaultUnsignedLimits, "MeV"),
            std::tuple("pi1_PE", defaultUnsignedLimits, "MeV"),
            std::tuple("pi1_PX", defaultSignedLimits, "MeV"),
            std::tuple("pi1_PY", defaultSignedLimits, "MeV"),
            std::tuple("pi1_PZ", defaultSignedLimits, "MeV"),
            std::tuple("pi1_M", defaultUnsignedLimits, "MeV"),
            std::tuple("pi1_ID", defaultSignedLimits, ""),

            std::tuple("pi2_IP_OWNPV", defaultUnsignedLimits, ""),
            std::tuple("pi2_IPCHI2_OWNPV", defaultUnsignedLimits, ""),
            std::tuple("pi2_P", defaultUnsignedLimits, "MeV"),
            std::tuple("pi2_PT", defaultUnsignedLimits, "MeV"),
            std::tuple("pi2_PE", defaultUnsignedLimits, "MeV"),
            std::tuple("pi2_PX", defaultSignedLimits, "MeV"),
            std::tuple("pi2_PY", defaultSignedLimits, "MeV"),
            std::tuple("pi2_PZ", defaultSignedLimits, "MeV"),
            std::tuple("pi2_M", defaultUnsignedLimits, "MeV"),
            std::tuple("pi2_ID", defaultSignedLimits, ""),

            std::tuple("pi3_IP_OWNPV", defaultUnsignedLimits, ""),
            std::tuple("pi3_IPCHI2_OWNPV", defaultUnsignedLimits, ""),
            std::tuple("pi3_P", defaultUnsignedLimits, "MeV"),
            std::tuple("pi3_PT", defaultUnsignedLimits, "MeV"),
            std::tuple("pi3_PE", defaultUnsignedLimits, "MeV"),
            std::tuple("pi3_PX", defaultSignedLimits, "MeV"),
            std::tuple("pi3_PY", defaultSignedLimits, "MeV"),
            std::tuple("pi3_PZ", defaultSignedLimits, "MeV"),
            std::tuple("pi3_M", defaultUnsignedLimits, "MeV"),
            std::tuple("pi3_ID", defaultSignedLimits, ""),

            std::tuple("pisoft_IP_OWNPV", defaultUnsignedLimits, ""),
            std::tuple("pisoft_IPCHI2_OWNPV", defaultUnsignedLimits, ""),
            std::tuple("pisoft_P", defaultUnsignedLimits, "MeV"),
            std::tuple("pisoft_PT", defaultUnsignedLimits, "MeV"),
            std::tuple("pisoft_PE", defaultUnsignedLimits, "MeV"),
            std::tuple("pisoft_PX", defaultSignedLimits, "MeV"),
            std::tuple("pisoft_PY", defaultSignedLimits, "MeV"),
            std::tuple("pisoft_PZ", defaultSignedLimits, "MeV"),
            std::tuple("pisoft_M", defaultUnsignedLimits, "MeV"),
            std::tuple("pisoft_ID", defaultSignedLimits, ""),

            std::tuple("DELTA_M", std::pair(139.3, 168.132), "MeV")}; // hard coded ew
}

void createWeightedRootFile(const std::string&              inFile,
                            const std::string&              treeName,
                            RooAbsPdf&                      signalModel,
                            RooAbsPdf&                      backgroundModel,
                            const std::vector<std::string>& fixedParams,
                            const std::string&              outTreePath)
{
    // Create a RooWorkspace for doing the sWeighting
    RooWorkspace workspace("sWeighting workspace");

    // Read in the data we want from the tree; first select which branches we want then read them into the workspace
    std::vector<std::string>               branches{};
    std::vector<std::pair<double, double>> ranges{};
    std::vector<std::string>               units{};
    for (const auto& c : desiredBranches()) {
        branches.push_back(std::get<0>(c));
        ranges.push_back(std::get<1>(c));
        units.push_back(std::get<2>(c));
    }
    readData(workspace, inFile, treeName, branches, ranges, units);

    // Create a combined model and perform sWeighting
    TTree* weightedData = sWeightData(workspace, signalModel, backgroundModel, fixedParams);

    // TODO Show the plot maybe

    // Write the sWeighted data to a new nTuple
    weightedData->SaveAs(outTreePath.c_str());
    delete weightedData;
}

} // namespace sWeighting
