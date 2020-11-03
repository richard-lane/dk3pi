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
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLine.h>

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
 * Create a plot of the fit
 */
static void massFitPlot(const RooDataSet&  dataset,
                        const RooArgSet&   params,
                        const std::string& observable,
                        RooAbsPdf&         combinedModel,
                        const RooAbsPdf&   signalModel,
                        const RooAbsPdf&   backgroundModel,
                        const char*        plotPath)
{
    // Create a canvas with appropriate pads for a main mass fit plot and a subsidiary residual plot
    std::unique_ptr<TCanvas> c{std::make_unique<TCanvas>("data", "data")};
    std::unique_ptr<TPad> massFitPad{std::make_unique<TPad>("mass fit", "mass fit", 0, 0.2, 1, 1)};
    std::unique_ptr<TPad> residualPad{std::make_unique<TPad>("residual", "pull", 0, 0, 1, 0.2)};
    massFitPad->SetBottomMargin(0.00001);
    massFitPad->SetBorderMode(0);
    residualPad->SetTopMargin(0.00001);
    residualPad->SetBottomMargin(0.35);
    residualPad->SetBorderMode(0);
    massFitPad->Draw();
    residualPad->Draw();

    // Create a RooPlot for our observable
    RooRealVar*              observableData = dynamic_cast<RooRealVar*>(&params[TString(observable)]);
    const int                numBins{70};
    std::unique_ptr<RooPlot> fitPlot(observableData->frame(RooFit::Title("Mass Fit"), RooFit::Bins(numBins)));

    // Add the model and data to the RooPlot
    // Add the combined model last, as this is the one that RooPlot will "remember"
    // i.e. it will use these for calculating pulls, residuals etc.
    dataset.plotOn(fitPlot.get());
    combinedModel.plotOn(
        fitPlot.get(), RooFit::Components(signalModel), RooFit::LineStyle(kDashed), RooFit::LineColor(kBlue));
    combinedModel.plotOn(
        fitPlot.get(), RooFit::Components(backgroundModel), RooFit::LineStyle(kDotted), RooFit::LineColor(kRed));
    combinedModel.plotOn(fitPlot.get(), RooFit::LineColor(kGreen));

    // Create a pull and a new frame to contain it
    RooHist*                 pullHist{fitPlot->pullHist()};
    std::unique_ptr<RooPlot> pullPlot(observableData->frame(RooFit::Title(" "), RooFit::Bins(numBins))); // Blank title
    pullPlot->addPlotable(pullHist, "P");
    pullPlot->GetYaxis()->SetNdivisions(10);

    // Increase the size of the x axis label and tick size
    pullPlot->GetXaxis()->SetTitleSize(0.15);
    pullPlot->GetXaxis()->SetLabelSize(0.15);

    // Draw the mass fit
    massFitPad->cd();
    fitPlot->Draw();

    // Draw the pull and a horizontal line
    TLine line(observableData->getMin(), 0, observableData->getMax(), 0);
    residualPad->cd();
    pullPlot->Draw();
    line.Draw();

    c->SaveAs(plotPath);
};

/*
 * Perform sWeighting
 *
 * The sPlot method requires us to fix all parameters of the model that are not yields once we have performed the
 * fit; these should be passed by name via fixedParameters
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
                                          const char*                     massFitPlotPath = nullptr,
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
    if (massFitPlotPath) {
        massFitPlot(data, *params, observable, combinedModel, signalModel, backgroundModel, massFitPlotPath);
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

std::unique_ptr<TTree> createSWeightedTree(const std::string&              inFile,
                                           const std::string&              treeName,
                                           std::vector<RootBranch>         branchInfo,
                                           RooAbsPdf&                      signalModel,
                                           RooAbsPdf&                      backgroundModel,
                                           const std::string&              observable,
                                           const std::vector<std::string>& fixedParams,
                                           const char*                     massFitPlot,
                                           const char*                     graphVizDiagram)
{

    // Read in the data we want from the tree; first select which branches we want then read the data from them into
    // a RooDataSet
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
