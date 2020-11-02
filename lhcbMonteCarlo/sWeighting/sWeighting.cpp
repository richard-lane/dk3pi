/*
 * Based on Jenny Smallwood's (Oxford) code that she sent as an example of D->K3Pi sWeighting
 */
#include "sWeighting.h"

#include <memory>
#include <tuple>

#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooStats/SPlot.h>
#include <TCanvas.h>
#include <TFile.h>

namespace sWeighting
{

/*
 * Read data from a ROOT file into a RooDataSet
 *
 * Need to pass in vectors of branch names, ranges (i.e. min/max allowed values) and units (e.g. "MeV^2")
 */
static RooDataSet readData(const std::string&                            rootFilePath,
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
    return RooDataSet(
        "data", treeName.c_str(), *branches, RooFit::ImportFromFile(rootFilePath.c_str(), treeName.c_str()));
}

/*
 * Perform sWeighting
 *
 * The sPlot method requires us to fix all parameters of the model that are not yields once we have performed the fit;
 * these should be passed by name via fixedParameters
 *
 * NB: doesn't check that sufficient/sensible parameters are fixed, since that's annoyingly hard
 *
 * if a mass fit plot C-string is provided then a plot of the mass fit will be created
 * if a graphViz diagram C-string is provided then a graph showing the model structure will be created
 *
 */
static std::unique_ptr<TTree> sWeightData(RooDataSet&                     data,
                                          RooAbsPdf&                      signalModel,
                                          RooAbsPdf&                      backgroundModel,
                                          const std::string&              observable,
                                          const std::vector<std::string>& fixedParameters,
                                          const char*                     massFitPlot     = nullptr,
                                          const char*                     graphVizDiagram = nullptr)
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

    // Create a graphviz .dot file to visualise the fitter, which might be useful
    if (graphVizDiagram) {
        combinedModel.graphVizTree(graphVizDiagram);
    }

    // Fit the model to data, using an extended likelihood fit
    combinedModel.fitTo(data, RooFit::Extended());

    // Create an argument set the contains our models' parameters
    // This allocates memory!
    RooArgSet* params{combinedModel.getVariables()};

    // Create a mass fit plot if we need to
    if (massFitPlot) {
        RooRealVar* observableData = dynamic_cast<RooRealVar*>(&(*params)[TString(observable)]);

        TCanvas* c     = new TCanvas("data", "data");
        RooPlot* frame = new RooPlot(
            "Mass Fit", "Mass Fit", *observableData, observableData->getMin(), observableData->getMax(), 100);
        data.plotOn(frame);
        combinedModel.plotOn(frame);
        combinedModel.plotOn(frame, RooFit::Components(signalModel), RooFit::LineStyle(kDashed));
        combinedModel.plotOn(frame, RooFit::Components(backgroundModel), RooFit::LineStyle(kDotted));
        frame->Draw();
        c->SaveAs(massFitPlot);
        delete frame;
        delete c;
    }

    // Fix the parameters that we were meant to fix
    for (const auto& paramName : fixedParameters) {
        dynamic_cast<RooRealVar*>(&(*params)[paramName.c_str()])->setConstant();
    }

    // Perform sPlot fit to find the number of signal and background events
    RooStats::SPlot("sData", "An sPlot", data, &combinedModel, RooArgList(numSignalEvents, numBackgroundEvents));

    delete params;
    return std::unique_ptr<TTree>(data.GetClonedTree());
}

std::unique_ptr<TTree> createWeightedRootFile(const std::string&              inFile,
                                              const std::string&              treeName,
                                              std::vector<RootBranch>         branchInfo,
                                              RooAbsPdf&                      signalModel,
                                              RooAbsPdf&                      backgroundModel,
                                              const std::string&              observable,
                                              const std::vector<std::string>& fixedParams,
                                              const char*                     massFitPlot,
                                              const char*                     graphVizDiagram)
{

    // Read in the data we want from the tree; first select which branches we want then read the data from them into a
    // RooDataSet
    std::vector<std::string>               branches{};
    std::vector<std::pair<double, double>> ranges{};
    std::vector<std::string>               units{};
    for (const auto& c : branchInfo) {
        branches.push_back(c.name);
        ranges.push_back(std::make_pair(c.min, c.max));
        units.push_back(c.units);
    }
    RooDataSet data = readData(inFile, treeName, branches, ranges, units);

    // Create a combined model and perform sWeighting
    return sWeightData(data, signalModel, backgroundModel, observable, fixedParams, massFitPlot, graphVizDiagram);
}

} // namespace sWeighting
