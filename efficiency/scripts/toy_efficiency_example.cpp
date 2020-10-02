/*
 * Read in some data from AmpGen, apply a toy efficiency function and plot Dalitz plots and 1d projections of invariant
 * masses or something
 */
#include <iostream>
#include <memory>
#include <random>
#include <string>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLegend.h>

#include "ReadRoot.h"
#include "bdt_reweighting.h"
#include "efficiency.h"
#include "efficiencyUtil.h"
#include "scriptUtils.h"
#include "util.h"

static std::vector<PhspPoint> parametrisePoints(const std::vector<dDecay_t>& decay)
{
    std::vector<PhspPoint> parametrisedPoints(decay.size());
    for (size_t i = 0; i < decay.size(); ++i) {
        parametrisedPoints[i] = parametrisation(decay[i]);
    }
    return parametrisedPoints;
}

static void plotHists(TH1D& truthHist, TH1D& mcHist, TH1D& detectedHist, TH1D& correctedHist, const char* path)
{
    std::unique_ptr<TCanvas> canvasPtr = std::make_unique<TCanvas>();
    canvasPtr->cd();
    canvasPtr->SetLeftMargin(0.15);

    correctedHist.SetLineColor(kBlue);
    mcHist.SetLineColor(kBlack);
    truthHist.SetLineColor(kGreen);
    detectedHist.SetLineColor(kRed);

    correctedHist.SetStats(false);
    mcHist.SetStats(false);
    truthHist.SetStats(false);
    detectedHist.SetStats(false);

    // Legend
    std::unique_ptr<TLegend> legendPtr = std::make_unique<TLegend>(0.15, 0.25);
    legendPtr->SetTextSize(0.03);
    legendPtr->AddEntry(&correctedHist, "Corrected", "l");
    legendPtr->AddEntry(&mcHist, "MC", "l");
    legendPtr->AddEntry(&truthHist, "Truth", "l");
    legendPtr->AddEntry(&detectedHist, "Detected", "l");

    correctedHist.Draw("SAME");
    mcHist.Draw("SAME");
    truthHist.Draw("SAME");
    detectedHist.Draw("SAME");
    legendPtr->Draw();
    canvasPtr->SaveAs(path);

    // If you want to visualise the actual efficiency
    TH1D* tmp = (TH1D*)detectedHist.Clone();
    tmp->SetTitle((std::string("Efficiency") + path).c_str());
    tmp->Divide(&truthHist);
    util::saveObjectToFile(tmp, (std::string("efficiency_") + path).c_str());
}

/*
 * Take a collection of truth phase space points, detected phsp points and weights + make plots of truth parameters,
 * detected parameters + weighted detected parameters
 */
static void makePlots(const std::vector<PhspPoint>& truth,
                      const std::vector<PhspPoint>& mc,
                      const std::vector<PhspPoint>& detected,
                      const std::vector<double>     weights,
                      const PhspBins&               binLimits,
                      const std::string&            title)
{
    size_t dimensionality = truth[0].size();
    assert(dimensionality == detected[0].size());
    assert(weights.size() == detected.size());

    for (size_t i = 0; i < dimensionality; ++i) {
        // Translate an index into a parametrisation label
        std::string paramLabel = "";
        switch (i) {
        case 0: paramLabel = "m(K- pi+)"; break;
        case 1: paramLabel = "m(pi+ pi+)"; break;
        case 2: paramLabel = "m(pi+ pi-)"; break;
        case 3: paramLabel = "m(K- pi+ pi-)"; break;
        case 4: paramLabel = "m(pi+ pi+ pi-)"; break;
        default:
            // We're in 5d
            std::cerr << "Expected 5d phsp, got dimension index" << i << std::endl;
            assert(false);
        }

        // Assign memory to histograms
        std::unique_ptr<TH1D> truthHist = std::make_unique<TH1D>((title + "truth_" + std::to_string(i)).c_str(),
                                                                 (title + "- truth " + std::to_string(i)).c_str(),
                                                                 binLimits[i].size() - 1,
                                                                 binLimits[i].data());
        std::unique_ptr<TH1D> mcHist    = std::make_unique<TH1D>((title + "mc_" + std::to_string(i)).c_str(),
                                                              (title + "- MC" + std::to_string(i)).c_str(),
                                                              binLimits[i].size() - 1,
                                                              binLimits[i].data());

        std::unique_ptr<TH1D> detectedHist = std::make_unique<TH1D>((title + "detected_" + std::to_string(i)).c_str(),
                                                                    (title + "- detected " + std::to_string(i)).c_str(),
                                                                    binLimits[i].size() - 1,
                                                                    binLimits[i].data());

        // It's this hist that determines that plot title
        std::unique_ptr<TH1D> correctedHist = std::make_unique<TH1D>(
            (title + "corrected_" + std::to_string(i)).c_str(),
            (title + " Reweighting Projection " + paramLabel + ";Inv Mass " + paramLabel + "/GeV; Count").c_str(),
            binLimits[i].size() - 1,
            binLimits[i].data());
        // Fill each histogram
        for (const PhspPoint& event : truth) {
            truthHist->Fill(event[i]);
        }
        for (const PhspPoint& event : mc) {
            mcHist->Fill(event[i]);
        }
        for (const PhspPoint& event : detected) {
            detectedHist->Fill(event[i]);
        }
        for (size_t j = 0; j < detected.size(); ++j) {
            if (std::isfinite(weights[j])) {
                correctedHist->Fill(detected[j][i], weights[j]);
            } else {
                std::cerr << "Warning: stupid weight " << weights[j] << " at index " << j << std::endl;
            }
        }
        plotHists(*truthHist, *mcHist, *detectedHist, *correctedHist, (title + std::to_string(i) + ".png").c_str());
    }
}

/*
 * Make an estimate to the efficiency using Chow Liu parametrisation of truth and mc, and reweight realData with these
 * efficiencies.
 *
 * Then plot
 */
static void chowLiu(const PhspBins&               bins,
                    const std::vector<PhspPoint>& truth,
                    const std::vector<PhspPoint>& mc,
                    const std::vector<PhspPoint>& realData)
{
    // Create the object used for making the efficiency parametrisation
    std::cout << "\n====Chow Liu ====" << std::endl;
    ChowLiuEfficiency EfficiencyCorrection(bins);

    // Add the events of both type
    std::cout << "Add truth events..." << std::flush;
    for (auto event : truth) {
        EfficiencyCorrection.addGeneratedEvent(event);
    }
    std::cout << "done" << std::endl;
    std::cout << "add detected events..." << std::flush;
    for (auto event : mc) {
        EfficiencyCorrection.addMCEvent(event);
    }
    std::cout << "done" << std::endl;

    // Perform efficiency parametrisation
    std::cout << "Perform efficiency parametrisation..." << std::flush;
    EfficiencyCorrection.efficiencyParametrisation();
    std::cout << "done" << std::endl;

    std::vector<double> chowLiuWeights(realData.size(), 0);
    for (size_t i = 0; i < chowLiuWeights.size(); ++i) {
        chowLiuWeights[i] = 1 / EfficiencyCorrection.value(realData[i]);
    }
    makePlots(truth, mc, realData, chowLiuWeights, bins, "chowLiu");
}

/*
 * Use BDT to make estimate of efficiencies using truth and MC, then plot realData with weights
 *
 * Uses bins for plotting
 */
static void bdt(const PhspBins&               bins,
                const std::vector<PhspPoint>& truth,
                const std::vector<PhspPoint>& mc,
                const std::vector<PhspPoint>& realData)
{
    // Take our data, train our BDT
    std::cout << "\n====BDT====" << std::endl;
    std::cout << "Training BDT..." << std::flush;
    PyObject* bdt = initBDT(truth, mc);
    std::cout << "done" << std::endl;

    std::cout << "Finding weights..." << std::flush;
    std::vector<double> weights = efficiency(bdt, realData, truth.size());
    std::cout << "done" << std::endl;

    makePlots(truth, mc, realData, weights, bins, "BDT");
}

/*
 * Read in a load of D->K3Pi events from a ROOT file (probably generated with AmpGen); this will serve as our
 * truth-level data
 *
 * Then prune these events to have "MC" and "real" datasets- these will follow the same distribution.
 *
 * Use MC/truth data to find an efficiency estimate, and plot truth-level and reweighted real data
 */
int main()
{
    // Read in data from a ROOT file that I generated with AmpGen to get out mock "truth-level" data
    std::string              ampgenRootFile("../../AmpGen/binning/Mixed.root");
    std::string              ampgenTreeName          = "DalitzEventList";
    std::vector<std::string> ampgenBranchNames       = {"_1_K~", "_2_pi#", "_3_pi#", "_4_pi~"};
    std::vector<std::string> ampgenMomentumPostfixes = {"_Px", "_Py", "_Pz", "_E"};
    std::unique_ptr<TFile>   tFile(new TFile(ampgenRootFile.c_str()));
    std::cout << "Read Ampgen data ...";
    ReadRoot RootData(tFile.get(), ampgenTreeName, ampgenBranchNames, ampgenMomentumPostfixes);
    std::cout << "done" << std::endl;

    // RNG that we need for applying efficiency to our data
    std::random_device rd;
    std::mt19937       generator(rd());

    // Function that we'll use as our efficiency
    const std::function<double(dDecay_t)> efficiencyFcn = awkwardEfficiency;

    std::cout << "Generate \"MC\" data.." << std::flush;
    std::vector<dDecay_t> trainingEvents = RootData.events;
    applyEfficiency(&generator, efficiencyFcn, trainingEvents);
    std::cout << "done" << std::endl;

    // Generate more MC data
    std::cout << "Generate \"real\" data..." << std::flush;
    std::vector<dDecay_t> detectedEvents = RootData.events;
    applyEfficiency(&generator, efficiencyFcn, detectedEvents);
    std::cout << "done" << std::endl;

    PhspBins Bins = findBins();

    // Find the phsp parametrisations of these datasets
    std::vector<PhspPoint> truth    = parametrisePoints(RootData.events);
    std::vector<PhspPoint> mc       = parametrisePoints(trainingEvents);
    std::vector<PhspPoint> realData = parametrisePoints(detectedEvents);

    // Make plots of truth/detected/reweighted events using each of our methods
    chowLiu(Bins, truth, mc, realData);
    bdt(Bins, truth, mc, realData);

    return 0;
}
