#include <memory>

#include "sWeighting.h"

#include <RooRealVar.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TTree.h>

// I stole these numbers from Jenny
// I don't know what i was thinking when i decided this was the thing to do
namespace JennysParams
{
static constexpr double lowDMass{139.3};
static constexpr double highDMass{168.0};

static constexpr double meanDMass{146.0};
static constexpr double meanDMassLoCutoff{144.0};
static constexpr double meanDMassHiCutoff{147.0};

static constexpr double dMassWidth{1.0};
static constexpr double dMassWidthLoCutoff{0.0001};
static constexpr double dMassWidthHiCutoff{5.0};

static constexpr double gamma{0.01};
static constexpr double gammaLoCutoff{-5.0};
static constexpr double gammaHiCutoff{5.0};

static constexpr double delta{1.0};
static constexpr double deltaLoCutoff{-10.0};
static constexpr double deltaHiCutoff{10.0};

static constexpr double a{0.01};
static constexpr double aLoCutoff{-10.0};
static constexpr double aHiCutoff{10.0};

static constexpr double b{0.01};
static constexpr double bLoCutoff{-10.0};
static constexpr double bHiCutoff{10.0};

static constexpr double c{2.0};
static constexpr double cLoCutoff{0.0001};
static constexpr double cHiCutoff{100.0};

static constexpr double massThreshhold{135};

// The fit observable: D*-D0 mass difference
static RooRealVar dMassVar("DELTA_M", "DELTA_M", JennysParams::lowDMass, JennysParams::highDMass, "MeV/c^2");

// Johnson signal fcn parameters
static RooRealVar meanDMassVar("Mean D Mass",
                               "Mean D Mass",
                               JennysParams::meanDMass,
                               JennysParams::meanDMassLoCutoff,
                               JennysParams::meanDMassHiCutoff);
static RooRealVar dMassWidthVar("D Mass Width",
                                "D Mass Width",
                                JennysParams::dMassWidth,
                                JennysParams::dMassWidthLoCutoff,
                                JennysParams::dMassWidthHiCutoff);
static RooRealVar
    gammaVar("Gamma", "Gamma", JennysParams::gamma, JennysParams::gammaLoCutoff, JennysParams::gammaHiCutoff);
static RooRealVar
    deltaVar("Delta", "Delta", JennysParams::delta, JennysParams::deltaLoCutoff, JennysParams::deltaHiCutoff);

// Dst-D0 background params
static RooRealVar deltaM0Var("Delta M0", "Delta M0", JennysParams::lowDMass);
static RooRealVar aVar("a", "a", JennysParams::a, JennysParams::aLoCutoff, JennysParams::aHiCutoff);
static RooRealVar bVar("b", "b", JennysParams::b, JennysParams::bLoCutoff, JennysParams::bHiCutoff);
static RooRealVar cVar("c", "c", JennysParams::c, JennysParams::cLoCutoff, JennysParams::cHiCutoff);

}; // namespace JennysParams

static inline RooJohnson jennysPromptSignalModel(void)
{

    return sWeighting::johnsonSignalModel(JennysParams::dMassVar,
                                          JennysParams::meanDMassVar,
                                          JennysParams::dMassWidthVar,
                                          JennysParams::gammaVar,
                                          JennysParams::deltaVar,
                                          JennysParams::massThreshhold);
}

static inline RooDstD0BG jennysPromptBackgroundModel(void)
{

    return sWeighting::promptBackgroundModel(
        JennysParams::dMassVar, JennysParams::deltaM0Var, JennysParams::aVar, JennysParams::bVar, JennysParams::cVar);
}

static void sPlotHist(const std::string &rootFile, const std::string &plotPath)
{
    // Read tree
    TFile *newFile = new TFile(rootFile.c_str());
    TTree *tree{nullptr};
    newFile->GetObject("RooTreeDataStore_workspace data_data;1", tree);
    assert(tree);

    // Set up buffers for reading data
    double signalWt{0};
    double backgroundWt{0};
    double deltaM{0};
    tree->SetBranchAddress("DELTA_M", &deltaM);
    tree->SetBranchAddress("numSignalEvents_sw", &signalWt);
    tree->SetBranchAddress("numBackgroundEvents_sw", &backgroundWt);

    // Create + fill hists
    TH1D *signal     = new TH1D("signal", "Reweighted Events;Delta_M /MeV; Weighted Count", 100, 135, 170);
    TH1D *background = new TH1D("bkg", "Reweighted Events;Delta_M /MeV; Weighted Count", 100, 135, 170);
    auto  numEvents{tree->GetEntries()};
    for (decltype(numEvents) i{0}; i < numEvents; ++i) {
        tree->GetEntry(i);
        signal->Fill(deltaM, signalWt);
        background->Fill(deltaM, backgroundWt);
    }

    signal->SetLineColor(kBlue);
    background->SetLineColor(kRed);

    TCanvas *c{new TCanvas()};
    c->cd();
    signal->SetStats(false);
    signal->Draw("HIST L SAME");
    background->Draw("HIST L SAME");

    // Legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(signal, "Signal");
    legend->AddEntry(background, "Background");
    legend->Draw("SAME");

    c->SaveAs(plotPath.c_str());

    delete c;
    delete legend;
    delete signal;
    delete background;
    delete newFile;
}

int main()
{
    // The file we want to read
    const std::string path{"./test_2011_RS_prompt.root"};
    const std::string treeName{"dk3pi/DecayTree"};

    // The parameters in our models that will be fixed after fitting
    const std::vector<std::string> paramsToFix{
        "a", "b", "c", "Gamma", "Delta", "Delta M0", "D Mass Width", "Mean D Mass"};

    // Create a signal model
    RooJohnson signalModel = jennysPromptSignalModel();

    // Create a background model
    RooDstD0BG backgroundModel = jennysPromptBackgroundModel();

    // Add a delta_m branch to our root file
    sWeighting::addDeltaMBranch(path, treeName, "D_M", "Dstar_M", "DELTA_M");

    // sWeight the data and output to a new ROOT file
    // The new DELTA_M branch gets written to a tree at /DecayTree/, not /dk3pi/DecayTree/ so i guess let's use
    std::string outFile{"newfile.root"};
    sWeighting::createWeightedRootFile(path, "DecayTree", signalModel, backgroundModel, paramsToFix, outFile);

    // Read in the new ROOT file and plot a histogram of reweighted mass difference
    sPlotHist(outFile, "reweighted.png");

    return 0;
}
