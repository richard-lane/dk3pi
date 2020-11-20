/*
 * Weight prompt real data to semileptonic
 */
#include "bdt_reweighting.h"
#include "scriptUtils.h"

#include <random>

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
 *
 * discard indices provided in the vector
 */
static std::vector<double> wts(const std::string&               rootFile,
                               const std::string&               treeName,
                               const std::string&               weightBranchName,
                               const std::vector<size_t>* const discard = nullptr)
{
    std::cout << "finding weights" << std::endl;
    // Create ROOT file and TTree
    TFile  f(rootFile.c_str());
    TTree* tree;
    f.GetObject(treeName.c_str(), tree);
    assert(tree);

    // Buffers
    double wt;
    tree->SetBranchAddress(weightBranchName.c_str(), &wt);

    // Init vector to the right length
    std::vector<double> weights{};

    // Iterate over tree filling in the right things
    // Unless we're meant to discard them
    size_t j{0};
    for (decltype(tree->GetEntries()) i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (!discard) {
            weights.push_back(wt);
        } else if ((size_t)i == (*discard)[j]) {
            j++;
        } else {
            weights.push_back(wt);
        }
    }

    delete tree;
    return weights;
}

/*
 * Take a D-> K pi1 pi2 pi3 root file, return a vector of phsppoints and a vector of the indices that i threw away
 */
static std::pair<std::vector<PhspPoint>, std::vector<size_t>>
phsp(const std::string& rootFile, const std::string& treeName, const bool prune = false)
{
    std::cout << "finding phsp" << std::endl;
    // Create ROOT file and TTree
    TFile  f(rootFile.c_str());
    TTree* tree;
    f.GetObject(treeName.c_str(), tree);
    assert(tree);

    // Buffers
    dDecay_t decay;
    tree->SetBranchAddress("D0_P0_PX", &decay.kParams.px);
    tree->SetBranchAddress("D0_P0_PY", &decay.kParams.py);
    tree->SetBranchAddress("D0_P0_PZ", &decay.kParams.pz);
    tree->SetBranchAddress("D0_P0_PE", &decay.kParams.energy);

    tree->SetBranchAddress("D0_P1_PX", &decay.pi1Params.px);
    tree->SetBranchAddress("D0_P1_PY", &decay.pi1Params.py);
    tree->SetBranchAddress("D0_P1_PZ", &decay.pi1Params.pz);
    tree->SetBranchAddress("D0_P1_PE", &decay.pi1Params.energy);

    tree->SetBranchAddress("D0_P2_PX", &decay.pi2Params.px);
    tree->SetBranchAddress("D0_P2_PY", &decay.pi2Params.py);
    tree->SetBranchAddress("D0_P2_PZ", &decay.pi2Params.pz);
    tree->SetBranchAddress("D0_P2_PE", &decay.pi2Params.energy);

    tree->SetBranchAddress("D0_P3_PX", &decay.pi3Params.px);
    tree->SetBranchAddress("D0_P3_PY", &decay.pi3Params.py);
    tree->SetBranchAddress("D0_P3_PZ", &decay.pi3Params.pz);
    tree->SetBranchAddress("D0_P3_PE", &decay.pi3Params.energy);

    // Init pair of empty vectors
    std::pair<std::vector<PhspPoint>, std::vector<size_t>> points{};

    // Iterate over tree filling in the right things
    for (decltype(tree->GetEntries()) i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);

        // Find out if this event is in the forbidden hole
        const double mass1 = invariantMass({decay.kParams, decay.pi1Params});
        const double mass2 = invariantMass({decay.pi1Params, decay.pi2Params});
        const double mass3 = invariantMass({decay.pi2Params, decay.pi3Params});
        const bool   eventInHole =
            (mass1 < 925.) && (mass1 > 875.) && (mass2 < 900.) && (mass2 > 700) && (mass3 < 900) && (mass3 > 700);

        // Don't add the point if we are pruning and the event is in the forbidden hole
        if (prune && eventInHole) {
            points.second.push_back(i);
        } else {
            points.first.push_back(parametrisation(decay));
        }
    }

    delete tree;
    return points;
}

static void
plotProjection(TH1D& promptHist, TH1D& slHist, TH1D& reweightedHist, const std::string& path, const std::string& xTitle)
{
    // Formatting
    promptHist.SetLineColor(kRed);
    slHist.SetLineColor(kGreen);
    reweightedHist.SetLineColor(kBlue);

    promptHist.SetStats(false);
    slHist.SetStats(false);
    reweightedHist.SetStats(false);

    // Scale hists
    promptHist.Scale(1 / promptHist.Integral());
    slHist.Scale(1 / slHist.Integral());
    reweightedHist.Scale(1 / reweightedHist.Integral());

    // Axis
    promptHist.GetYaxis()->SetTitle("Relative counts");
    promptHist.GetYaxis()->SetRangeUser(0., 0.045);
    slHist.GetYaxis()->SetRangeUser(0., 0.045);
    promptHist.GetXaxis()->SetTitle(xTitle.c_str());

    // Save the plots
    util::LegendParams_t legend{0.7, 0.9, 0.7, 0.9};
    util::saveObjectsToFile<TH1D>({&promptHist, &slHist, &reweightedHist},
                                  {"HIST E", "SAME HIST E", "SAME HIST E"},
                                  {"Prompt", "SL", "Reweighted Prompt"},
                                  path.c_str(),
                                  legend);
}

/*
 * Split a vector v into the first and second half
 *
 * Does something if the vector has an odd length
 */
template <typename T> static std::pair<std::vector<T>, std::vector<T>> splitVector(const std::vector<T>& v)
{
    const size_t halfSize = v.size() / 2;

    auto firstHalf  = std::vector<T>(v.begin(), v.begin() + halfSize);
    auto secondHalf = std::vector<T>(v.begin() + halfSize, v.end());

    return {firstHalf, secondHalf};
}

int main()
{
    // Prompt ROOT file (should contain a branch containing event weights)
    const std::string promptFile{"cut_wg_rs_prompt.root"};
    const std::string promptWeightFile("rs_weights.root");
    const std::string promptTree{"DecayTree"};
    const std::string promptWeightBranch{"numSignalEvents_sw"};

    // Semileptonic ROOT file (should contain a branch containing event weights)
    const std::string semileptonicFile{"cut_wg_rs_sl.root"};
    const std::string slWeightFile{"sl_weights.root"};
    const std::string semileptonicTree{"DecayTree"};
    const std::string semileptonicWeightBranch{"numSignalEvents_sw"};

    // Split the prompt data into two sets; this assumes they're distributed randomly in phase space
    auto promptPoints{phsp(promptFile, promptTree, true)};
    auto promptWeightArray{wts(promptWeightFile, promptTree, promptWeightBranch, &promptPoints.second)};

    auto promptDataHalves       = splitVector(promptPoints.first);
    auto firstHalfOfPromptData  = promptDataHalves.first;
    auto secondHalfOfPromptData = promptDataHalves.second;

    auto promptWeightHalves        = splitVector(promptWeightArray);
    auto firstHalfOfPromptWeights  = promptWeightHalves.first;
    auto secondHalfOfPromptWeights = promptWeightHalves.second;

    // Split the SL data into two sets; this assumes they're distributed randomly in phase space
    auto semileptonicPoints{phsp(semileptonicFile, semileptonicTree, true)};
    auto slWeightArray{wts(slWeightFile, semileptonicTree, semileptonicWeightBranch, &semileptonicPoints.second)};

    auto slDataHalves       = splitVector(semileptonicPoints.first);
    auto firstHalfOfSLData  = slDataHalves.first;
    auto secondHalfOfSLData = slDataHalves.second;

    auto slWeightHalves        = splitVector(slWeightArray);
    auto firstHalfOfSLWeights  = slWeightHalves.first;
    auto secondHalfOfSLWeights = slWeightHalves.second;

    // Train the BDT on the first half of the data
    std::cout << "Training BDT" << std::endl;
    PyObject* bdt = initBDT(firstHalfOfSLData, firstHalfOfPromptData, &firstHalfOfSLWeights, &firstHalfOfPromptWeights);

    // Reweight the second half of the prompt data to look like the prompt data
    std::cout << "Finding weights" << std::endl;
    auto efficiencyWeights{efficiency(bdt, secondHalfOfPromptData, semileptonicPoints.first.size())};

    // Calculate the weights that we need
    std::vector<double> prompt2SLweights = secondHalfOfPromptWeights;
    for (size_t i = 0; i < prompt2SLweights.size(); ++i) {
        prompt2SLweights[i] *= efficiencyWeights[i];
    }

    // Create histograms of semileptonic, prompt + the prompt half that has been reweighted
    constexpr int d{5};
    std::string   titles[d]{"proj0.png", "proj1.png", "proj2.png", "proj3.png", "proj4.png"};
    std::string   labels[5]{
        "M(K pi1) /MeV", "M(pi1 pi2) /MeV", "M(pi2 pi3) /MeV", "M(K pi1 pi2) /MeV", "M(pi1 pi2 pi3) /MeV"};
    double low[d]{400, 200, 200, 700, 400};
    double high[d]{1600, 1400, 1400, 1800, 1600};
    size_t nBins{100};

    for (int i = 0; i < d; ++i) {
        TH1D promptHist{plotProjection(
            secondHalfOfPromptData, secondHalfOfPromptWeights, i, "Phsp Projection", nBins, low[i], high[i])};
        TH1D semileptonicHist{
            plotProjection(secondHalfOfSLData, secondHalfOfSLWeights, i, "semileptonic", nBins, low[i], high[i])};
        TH1D reweightedPrompt{
            plotProjection(secondHalfOfPromptData, prompt2SLweights, i, "Reweighted", nBins, low[i], high[i])};
        plotProjection(promptHist, semileptonicHist, reweightedPrompt, titles[i], labels[i]);
    }

    // Create slices for something
    HistogramSlices promptSlices("Prompt Hist Slice;m(K\\pi_1);relative counts",
                                 10,
                                 100,
                                 std::make_pair(low[0], high[0]),
                                 std::make_pair(low[1], high[1]),
                                 0,
                                 1);
    promptSlices.add(secondHalfOfPromptData, &secondHalfOfPromptWeights);
    promptSlices.setColour(kRed);

    HistogramSlices slSlices("SL Hist Slice;m(K\\pi_1);relative counts",
                             10,
                             100,
                             std::make_pair(low[0], high[0]),
                             std::make_pair(low[1], high[1]),
                             0,
                             1);
    slSlices.add(semileptonicPoints.first, &slWeightArray);
    slSlices.setColour(kGreen);

    HistogramSlices reweightedSlices("Reweighted Hist Slice ;m(K\\pi_1);relative counts",
                                     10,
                                     100,
                                     std::make_pair(low[0], high[0]),
                                     std::make_pair(low[1], high[1]),
                                     0,
                                     1);
    reweightedSlices.add(secondHalfOfPromptData, &prompt2SLweights);
    reweightedSlices.setColour(kBlue);

    std::vector<HistogramSlices> slices{promptSlices, slSlices, reweightedSlices};
    plotSlices("slices",
               slices,
               {"HIST E SAME", "HIST E SAME", "HIST E SAME"},
               {"Prompt", "SL", "Reweighted"},
               {-0.0002, 0.015});
}
