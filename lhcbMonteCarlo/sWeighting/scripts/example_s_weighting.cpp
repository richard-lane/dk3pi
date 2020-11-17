#include <memory>

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
        new TH1D("signal", "Reweighted Events;M(D)/MeV; Weighted Count", 100, axisLimits.first, axisLimits.second);
    TH1D* background =
        new TH1D("bkg", "Reweighted Events;M(D)/MeV; Weighted Count", 100, axisLimits.first, axisLimits.second);
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
    const double loDeltaM{134.7};
    const double hiDeltaM{170.1};
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

    // sWeight the data
    // The new DELTA_M branch gets written to a tree at /DecayTree/, not /dk3pi/DecayTree/ so i guess let's use that
    sWeighting::Observable_t observable(observableName, loDeltaM, hiDeltaM, "MeV/c^2");
    sWeighting::findSWeights(inFile,
                             outFile,
                             "DecayTree",
                             signalModel,
                             backgroundModel,
                             500000,
                             500000,
                             observable,
                             paramsToFix,
                             outMassFitPlot.c_str(),
                             graph.c_str());

    // Read in the new ROOT file and plot a histogram of reweighted mass difference
    sPlotHist(outFile, outImg.c_str(), "DecayTree", "DELTA_M", {loDeltaM, hiDeltaM});
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

    // Create a signal model
    // These RooRealVars need to have the same scope as the models that use them, otherwise the models will segfault
    std::string  observableName{"D0_M"};
    const double lowDMass{1803.};
    const double hiDMass{1925.};
    RooRealVar   dMassVar(observableName.c_str(), observableName.c_str(), lowDMass, hiDMass, "MeV/c^2");
    RooRealVar   meanDMassVar("Mean D Mass", "Mean D Mass", 1865.0, 1850.0, 1880.0);
    RooRealVar   dMassWidthVar("D Mass Width", "D Mass Width", 10.0, 0.0001, 50.0);
    RooGaussian  signalModel("Signal Gaussian", "Signal PDF: Gaussian Function", dMassVar, meanDMassVar, dMassWidthVar);

    // Create a background model
    RooRealVar    aVar("a", "a", 0.005, 0.0, 0.01);
    RooRealVar    bVar("b", "b", 0.0, -10.0, 10.0);
    RooPolynomial backgroundModel("SL Background", "Background PDF: Polynomial", dMassVar, RooArgList(aVar, bVar));

    // Create a new ROOT file with the tree not in a directory
    // Hack for now until i get the sWeighting to properly deal with directories; TODO

    // The parameters in our models that will be fixed after fitting
    const std::vector<std::string> paramsToFix{"a", "b", observableName, "D Mass Width", "Mean D Mass"};

    // sWeight the data
    sWeighting::Observable_t observable(observableName, lowDMass, hiDMass, "MeV/c^2");
    sWeighting::findSWeights(inFile,
                             outFile,
                             "DecayTree",
                             signalModel,
                             backgroundModel,
                             500000,
                             500000,
                             observable,
                             paramsToFix,
                             outMassFitPlot.c_str(),
                             graph.c_str());

    // Read in the new ROOT file and plot a hist of reweighted D mass
    sPlotHist(outFile, outImg.c_str(), "DecayTree", observableName, {lowDMass, hiDMass});
}

int main()
{
    const std::string rsPath{"wg_rs_prompt.root"};
    promptFit(rsPath, "rs_weights.root", "rs.png", "rsMassFit.png", "rsGraph.dot");

    // const std::string wsPath{"./test_2011_WS_prompt.root"};
    // promptFit(wsPath, "newWS.root", "ws.png", "wsMassFit.png", "wsGraph.dot");

    const std::string slPath{"wg_rs_sl.root"};
    semiLeptonicFit(slPath, "sl_weights.root", "sl.png", "slMassFit.png", "sl.dot");

    return 0;
}
