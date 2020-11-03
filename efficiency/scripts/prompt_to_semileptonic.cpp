/*
 * Weight prompt real data to semileptonic
 */
#include "bdt_reweighting.h"
#include "scriptUtils.h"

#include <TFile.h>
#include <TH1D.h>

/*
 * Create a hist for the projection of the index'th dimension of points
 */
static TH1D plotProjection(const std::vector<PhspPoint>& points,
                           const std::vector<double>&    weights,
                           const unsigned int            index,
                           const std::string&            title,
                           const size_t                  nBins,
                           const double                  low,
                           const double                  high)
{
    assert(index < points[0].size());
    assert(points.size() == weights.size());

    TH1D hist{TH1D(title.c_str(), title.c_str(), nBins, low, high)};

    for (size_t i = 0; i < points.size(); ++i) {
        hist.Fill(points[i][index], weights[i]);
    }

    return hist;
}

/*
 * Take a D-> K pi1 pi2 pi3 root file, return a vector of weights
 */
static std::vector<double>
wts(const std::string& rootFile, const std::string& treeName, const std::string& weightBranchName)
{
    // Create ROOT file and TTree
    TFile  f(rootFile.c_str());
    TTree* tree;
    f.GetObject(treeName.c_str(), tree);
    assert(tree);

    // Buffers
    double wt;
    tree->SetBranchAddress(weightBranchName.c_str(), &wt);

    // Init vector to the right length
    std::vector<double> weights(tree->GetEntries());

    // Iterate over tree filling in the right things
    for (decltype(tree->GetEntries()) i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        weights[i] = wt;
    }

    delete tree;
    return weights;
}

/*
 * Take a D-> K pi1 pi2 pi3 root file, return a vector of phsppoints
 */
static std::vector<PhspPoint> phsp(const std::string& rootFile, const std::string& treeName)
{
    // Create ROOT file and TTree
    TFile  f(rootFile.c_str());
    TTree* tree;
    f.GetObject(treeName.c_str(), tree);
    assert(tree);

    // Buffers
    dDecay_t decay;
    tree->SetBranchAddress("K_PX", &decay.kParams.px);
    tree->SetBranchAddress("K_PY", &decay.kParams.py);
    tree->SetBranchAddress("K_PZ", &decay.kParams.pz);
    tree->SetBranchAddress("K_PE", &decay.kParams.energy);

    tree->SetBranchAddress("pi1_PX", &decay.pi1Params.px);
    tree->SetBranchAddress("pi1_PY", &decay.pi1Params.py);
    tree->SetBranchAddress("pi1_PZ", &decay.pi1Params.pz);
    tree->SetBranchAddress("pi1_PE", &decay.pi1Params.energy);

    tree->SetBranchAddress("pi2_PX", &decay.pi2Params.px);
    tree->SetBranchAddress("pi2_PY", &decay.pi2Params.py);
    tree->SetBranchAddress("pi2_PZ", &decay.pi2Params.pz);
    tree->SetBranchAddress("pi2_PE", &decay.pi2Params.energy);

    tree->SetBranchAddress("pi3_PX", &decay.pi3Params.px);
    tree->SetBranchAddress("pi3_PY", &decay.pi3Params.py);
    tree->SetBranchAddress("pi3_PZ", &decay.pi3Params.pz);
    tree->SetBranchAddress("pi3_PE", &decay.pi3Params.energy);

    // Init vector to the right length
    std::vector<PhspPoint> points(tree->GetEntries());

    // Iterate over tree filling in the right things
    for (decltype(tree->GetEntries()) i = 0; i < tree->GetEntries(); ++i) {
        // Work out its phsp params, insert into the vector
        tree->GetEntry(i);
        points[i] = parametrisation(decay);
    }

    delete tree;
    return points;
}

static void plotProjection(TH1D& promptHist, TH1D& slHist, const std::string& path, const std::string& xTitle)
{
    // Formatting
    promptHist.SetLineColor(kRed);
    slHist.SetLineColor(kGreen);
    promptHist.SetStats(false);
    slHist.SetStats(false);

    // Scale hists
    promptHist.Scale(1 / promptHist.Integral());
    slHist.Scale(1 / slHist.Integral());

    // Axis
    promptHist.GetYaxis()->SetTitle("Relative counts");
    promptHist.GetXaxis()->SetTitle(xTitle.c_str());

    // Save the plots
    util::LegendParams_t legend{0.7, 0.9, 0.7, 0.9};
    util::saveObjectsToFile<TH1D>({&promptHist, &slHist}, {"HIST", "SAME HIST"}, {"Prompt", "SL"}, path.c_str(), legend);
}

int main()
{
    // Prompt ROOT file (should contain a branch containing event weights)
    const std::string promptFile{"newRS.root"};
    const std::string promptTree{"RooTreeDataStore_data_DecayTree"};
    const std::string promptWeightBranch{"numSignalEvents_sw"};

    // Semileptonic ROOT file (should contain a branch containing event weights)
    const std::string semileptonicFile{"newSL.root"};
    const std::string semileptonicTree{"RooTreeDataStore_data_DecayTree"};
    const std::string semileptonicWeightBranch{"numSignalEvents_sw"};

    // Create vectors of phase space points
    auto promptPoints{phsp(promptFile, promptTree)};
    auto semileptonicPoints{phsp(semileptonicFile, semileptonicTree)};

    // Create vectors of weights
    auto promptWeights{wts(promptFile, promptTree, promptWeightBranch)};
    auto semileptonicWeights{wts(semileptonicFile, semileptonicTree, semileptonicWeightBranch)};

    // Reweight prompt to semileptonic

    // Create histograms
    const size_t nBins{100};
    const double low{200};
    const double high{1400};

    TH1D promptHist{plotProjection(promptPoints, promptWeights, 1, "Phsp Projection", nBins, low, high)};
    TH1D semileptonicHist{plotProjection(semileptonicPoints, semileptonicWeights, 1, "semileptonic", nBins, low, high)};

    plotProjection(promptHist, semileptonicHist, "proj1.png", "M(pi1 pi2) /MeV");
}
