/*
 * Based on Jenny Smallwood's (Oxford) code that she sent as an example of D->K3Pi sWeighting
 */
#include "sWeighting.h"

#include <memory>
#include <tuple>

#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooArgSet.h>
#include <RooDataHist.h>
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
 * Read observable data from a ROOT file into a RooDataSet
 *
 */
static RooDataSet readData(const std::string& rootFilePath, const std::string& treeName, const Observable_t& observable)
{
    std::cout << "Reading data from" << rootFilePath << std::endl;
    RooArgSet branch{};
    assert(branch.addClone(RooRealVar(observable.name.c_str(),
                                      observable.name.c_str(),
                                      observable.range.first,
                                      observable.range.second,
                                      observable.units.c_str())));

    // Create a RooDataSet, telling it compress data into a Tree
    RooAbsData::setDefaultStorageType(RooAbsData::Tree);
    return RooDataSet("data", treeName.c_str(), branch, RooFit::ImportFromFile(rootFilePath.c_str(), treeName.c_str()));
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
    std::unique_ptr<TPad>    massFitPad{std::make_unique<TPad>("mass fit", "mass fit", 0, 0.2, 1, 1)};
    std::unique_ptr<TPad>    residualPad{std::make_unique<TPad>("residual", "pull", 0, 0, 1, 0.2)};
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
    pullPlot->GetYaxis()->SetRangeUser(-30, 30);

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
 * Returns the params after fitting
 */
static std::unique_ptr<RooArgSet>
massFit(RooDataSet& data, RooAddPdf& combinedModel, const int numCPU, const std::string& observable)
{
    // Construct a dataset containing only our observable- this is all we need to perform the fit, and RooDataHist will
    // crash if passed a dataset containing too many vars
    const RooArgSet             fullArgSet = *data.get();
    RooArgSet                   observableArgSet(fullArgSet[TString(observable)]);
    std::unique_ptr<RooDataSet> minimalData(
        dynamic_cast<RooDataSet*>(data.reduce(RooFit::SelectVars(observableArgSet))));

    // Bin data TODO
    std::cout << "binning data" << std::endl;
    std::unique_ptr<RooDataHist> hist(minimalData->binnedClone());

    // Fit the model to data, using a multithreaded extended likelihood fit
    combinedModel.fitTo(*hist, RooFit::Extended(), RooFit::NumCPU(numCPU));

    // Create an argument set the contains our models' parameters
    std::unique_ptr<RooArgSet> params{combinedModel.getVariables()};

    return params;
}

/*
 * Perform sWeighting
 *
 * The sPlot method requires us to fix all parameters of the model that are not yields once we have performed the
 * fit; these should be passed by name via fixedParameters
 *
 * NB: doesn't check that sufficient/sensible parameters are fixed, since that's annoyingly hard
 *
 * if a mass fit plot C-string is provided then a plot of the mass fit will be created
 *
 * if a graphViz diagram C-string is provided then a graph showing the model structure will be created
 *
 * Returns a TTree containing the observable + some sWeighting branches
 *
 */
static std::unique_ptr<TTree> sWeightData(RooDataSet&                     data,
                                          RooAbsPdf&                      signalModel,
                                          RooAbsPdf&                      backgroundModel,
                                          const int                       expectedNumSignal,
                                          const int                       expectedNumBackground,
                                          const std::string&              observable,
                                          const std::vector<std::string>& fixedParameters,
                                          const int                       numCPU,
                                          const char*                     massFitPlotPath = nullptr,
                                          const char*                     graphVizDiagram = nullptr)
{
    // Create a combined model
    RooRealVar numSignalEvents(
        "numSignalEvents", "Signal Events", expectedNumSignal, 0, expectedNumSignal + expectedNumBackground);
    RooRealVar numBackgroundEvents("numBackgroundEvents",
                                   "Background Events",
                                   expectedNumBackground,
                                   0,
                                   expectedNumSignal + expectedNumBackground);
    RooAddPdf  combinedModel("combined model",
                            "Combined Signal + Background Models",
                            RooArgList(signalModel, backgroundModel),
                            RooArgList(numSignalEvents, numBackgroundEvents));

    // Create a graphviz .dot file to visualise the fitter, which might be useful
    if (graphVizDiagram) {
        combinedModel.graphVizTree(graphVizDiagram);
    }

    // Perform mass fit + return the parameters
    std::cout << "Performing mass fit" << std::endl;
    std::unique_ptr<RooArgSet> params = massFit(data, combinedModel, numCPU, observable);

    // Create a mass fit plot if we need to
    if (massFitPlotPath) {
        massFitPlot(data, *params, observable, combinedModel, signalModel, backgroundModel, massFitPlotPath);
    }

    // Fix the parameters that we were meant to fix before sWeighting
    for (const auto& paramName : fixedParameters) {
        dynamic_cast<RooRealVar*>(&(*params)[paramName.c_str()])->setConstant();
    }

    // Perform sPlot fit to find the number of signal and background events
    std::cout << "Creating sPlot" << std::endl;
    RooStats::SPlot("sData", "An sPlot", data, &combinedModel, RooArgList(numSignalEvents, numBackgroundEvents));

    return std::unique_ptr<TTree>(data.GetClonedTree());
}

void findSWeights(const std::string&              inFile,
                  const std::string&              outFile,
                  const std::string&              treeName,
                  RooAbsPdf&                      signalModel,
                  RooAbsPdf&                      backgroundModel,
                  const int                       expectedNumSignal,
                  const int                       expectedNumBackground,
                  const Observable_t&             observable,
                  const std::vector<std::string>& fixedParams,
                  const char*                     massFitPlot,
                  const char*                     graphVizDiagram,
                  const int                       numCPU)
{
    // Read data from the ROOT file into a RooDataSet
    RooDataSet data = readData(inFile, treeName, observable);

    // Create a combined model and perform sWeighting
    std::unique_ptr<TTree> sWeightingTree = sWeightData(data,
                                                        signalModel,
                                                        backgroundModel,
                                                        expectedNumSignal,
                                                        expectedNumBackground,
                                                        observable.name,
                                                        fixedParams,
                                                        numCPU,
                                                        massFitPlot,
                                                        graphVizDiagram);

    // Write to file
    sWeightingTree->SetName(treeName.c_str());
    sWeightingTree->SaveAs(outFile.c_str(), "RECREATE");
}

} // namespace sWeighting
