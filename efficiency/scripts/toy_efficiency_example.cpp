/*
 * Read in some data from AmpGen, apply a toy efficiency function and plot Dalitz plots and 1d projections of invariant
 * masses or something
 */
#include <memory>
#include <random>
#include <string>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLegend.h>

#include "ReadRoot.h"
#include "efficiency.h"
#include "efficiencyUtil.h"
#include "toyStudy.h"
#include "util.h"

PhspPoint parametrisation(const dDecay_t& decay)
{
    // 5d phsp
    PhspPoint point = PhspPoint(5);

    // Use invariant masses m12, m23, m34, m123, m234
    point[0] = invariantMass(std::vector<kinematicParams_t>{decay.kParams, decay.pi1Params});
    point[1] = invariantMass(std::vector<kinematicParams_t>{decay.pi1Params, decay.pi2Params});
    point[2] = invariantMass(std::vector<kinematicParams_t>{decay.pi2Params, decay.pi3Params});
    point[3] = invariantMass(std::vector<kinematicParams_t>{decay.kParams, decay.pi1Params, decay.pi2Params});
    point[4] = invariantMass(std::vector<kinematicParams_t>{decay.pi1Params, decay.pi2Params, decay.pi3Params});

    return point;
}

/*
 * Plot 1d projections of true, detected + reconstructed data
 */
void correctionPlot(const size_t                 param,
                    const std::vector<double>&   binLimits,
                    const std::vector<dDecay_t>& trueEvents,
                    const std::vector<dDecay_t>& detectedEvents,
                    const ChowLiuEfficiency&     EfficiencyCorrection,
                    const std::string&           path)
{
    // On the same canvas, draw histograms of the truth events, the measured events + the reweighted measured events
    std::string              title     = "Example Projection: " + path;
    std::unique_ptr<TCanvas> canvasPtr = std::make_unique<TCanvas>();
    canvasPtr->cd();
    canvasPtr->SetLeftMargin(0.15);

    std::unique_ptr<TH1D> truth =
        std::make_unique<TH1D>("truth", title.c_str(), binLimits.size() - 1, binLimits.data());
    std::unique_ptr<TH1D> detected =
        std::make_unique<TH1D>("detected", "detected", binLimits.size() - 1, binLimits.data());
    std::unique_ptr<TH1D> corrected =
        std::make_unique<TH1D>("corrected", "corrected", binLimits.size() - 1, binLimits.data());

    // Truth histogram is easy to fill
    std::cout << "Fill truth histogram..." << std::flush;
    for (auto truthEvent : trueEvents) {
        truth->Fill(parametrisation(truthEvent)[param]);
    }
    std::cout << "done" << std::endl;

    // Detected histogram is also easy to fill
    // Corrected histogram is filled with weights 1/efficiency for each point
    std::cout << "Calculate + fill corrected + detected histograms..." << std::flush;
    for (auto detectedEvent : detectedEvents) {
        detected->Fill(parametrisation(detectedEvent)[param]);
        double weight = 1 / EfficiencyCorrection.value(parametrisation(detectedEvent));
        corrected->Fill(parametrisation(detectedEvent)[param], weight);
    }
    std::cout << "done" << std::endl;

    truth->SetLineColor(kGreen);
    corrected->SetLineColor(kBlue);
    detected->SetLineColor(kRed);

    // Legend
    std::unique_ptr<TLegend> legendPtr = std::make_unique<TLegend>(0.15, 0.25);
    legendPtr->SetTextSize(0.03);

    legendPtr->AddEntry(truth.get(), "Truth", "l");
    truth->Draw("SAME");

    legendPtr->AddEntry(corrected.get(), "Corrected", "l");
    corrected->Draw("SAME");

    legendPtr->AddEntry(detected.get(), "Detected", "l");
    detected->Draw("SAME");

    legendPtr->Draw();
    canvasPtr->SaveAs(path.c_str());

    // If you want to visualise the actual efficiency
    TH1D tmp = *detected;
    tmp.SetTitle(("Efficiency" + path).c_str());
    tmp.Divide(truth.get());
    util::saveObjectToFile(&tmp, ("efficiency_" + path).c_str());
}

void dalitzPlot(const std::vector<dDecay_t>&     events,
                const std::pair<double, double>& axisLimits,
                const std::string&               title,
                const std::string&               path)
{
    size_t              numEvents = events.size();
    std::vector<double> m12(numEvents); // m(Kpi1)
    std::vector<double> m13(numEvents); // m(Kpi2)
    for (size_t i = 0; i < numEvents; ++i) {
        m12[i] = invariantMass(std::vector<kinematicParams_t>{events[i].kParams, events[i].pi1Params});
        m13[i] = invariantMass(std::vector<kinematicParams_t>{events[i].kParams, events[i].pi2Params});
    }

    // Plot them on a stupid, annoying, impossible ROOT histogram
    size_t                numBins = 100;
    double                min     = axisLimits.first;
    double                max     = axisLimits.second;
    std::unique_ptr<TH2D> m12m13Hist(
        new TH2D("dalitz", std::string(title + ";m12;m13").c_str(), numBins, min, max, numBins, min, max));
    m12m13Hist->FillN(numEvents, m12.data(), m13.data(), nullptr);

    //  Save Dalitz plots to file
    util::saveObjectToFile(m12m13Hist.get(), path.c_str(), "COLZ");
}

void histograms(const std::vector<dDecay_t>&     allEvents,
                const std::vector<dDecay_t>&     detectedEvents,
                const std::pair<double, double>& axisLimits,
                const std::string&               title,
                const std::string&               path)
{
    size_t              nAllEvents = allEvents.size();
    std::vector<double> allM12(nAllEvents); // m(Kpi1)
    std::vector<double> allM13(nAllEvents); // m(Kpi2)
    for (size_t i = 0; i < nAllEvents; ++i) {
        allM12[i] = invariantMass(std::vector<kinematicParams_t>{allEvents[i].kParams, allEvents[i].pi1Params});
        allM13[i] = invariantMass(std::vector<kinematicParams_t>{allEvents[i].kParams, allEvents[i].pi2Params});
    }

    size_t              nDetEvents = detectedEvents.size();
    std::vector<double> detM12(nAllEvents); // m(Kpi1)
    std::vector<double> detM13(nAllEvents); // m(Kpi2)
    for (size_t i = 0; i < nAllEvents; ++i) {
        detM12[i] =
            invariantMass(std::vector<kinematicParams_t>{detectedEvents[i].kParams, detectedEvents[i].pi1Params});
        detM13[i] =
            invariantMass(std::vector<kinematicParams_t>{detectedEvents[i].kParams, detectedEvents[i].pi2Params});
    }

    size_t                numBins = 100;
    std::unique_ptr<TH1D> allM12Hist(
        new TH1D("m12", (title + " m12").c_str(), numBins, axisLimits.first, axisLimits.second));
    std::unique_ptr<TH1D> allM13Hist(
        new TH1D("m23", (title + " m13").c_str(), numBins, axisLimits.first, axisLimits.second));

    std::unique_ptr<TH1D> detM12Hist(
        new TH1D("detm12", (title + " m12").c_str(), numBins, axisLimits.first, axisLimits.second));
    std::unique_ptr<TH1D> detM13Hist(
        new TH1D("detm23", (title + " m13").c_str(), numBins, axisLimits.first, axisLimits.second));

    allM12Hist->FillN(nAllEvents, allM12.data(), nullptr);
    allM13Hist->FillN(nAllEvents, allM13.data(), nullptr);
    detM12Hist->FillN(nDetEvents, detM12.data(), nullptr);
    detM13Hist->FillN(nDetEvents, detM13.data(), nullptr);

    std::unique_ptr<TCanvas> c(new TCanvas());
    c->cd();
    allM12Hist->Draw();
    detM12Hist->Draw("SAME");
    c->SaveAs((std::string("m12_") + path).c_str());

    std::unique_ptr<TCanvas> c2(new TCanvas());
    c2->cd();
    allM13Hist->Draw();
    detM13Hist->Draw("SAME");
    c2->SaveAs((std::string("m13_") + path).c_str());

    allM13Hist->Draw();
    detM13Hist->Draw("SAME");
}

int main()
{
    // Read in an AmpGen file
    // AmpGen
    std::string              ampgenRootFile("../../AmpGen/binning/Mixed.root");
    std::string              ampgenTreeName          = "DalitzEventList";
    std::vector<std::string> ampgenBranchNames       = {"_1_K~", "_2_pi#", "_3_pi#", "_4_pi~"};
    std::vector<std::string> ampgenMomentumPostfixes = {"_Px", "_Py", "_Pz", "_E"};
    std::unique_ptr<TFile>   tFile(new TFile(ampgenRootFile.c_str()));
    std::cout << "Read Ampgen data ...";
    ReadRoot RootData(tFile.get(), ampgenTreeName, ampgenBranchNames, ampgenMomentumPostfixes);
    std::cout << "done" << std::endl;

    // Create a lambda fcn that will reject all events with m12 below some value that i can read off the Dalitz plot
    // when it comes to it
    std::random_device rd;
    std::mt19937       generator(rd());

    // An efficiency with no correlations that gets recovered nicely
    auto simpleEfficiency = [](const dDecay_t& event) {
        (void)event;
        return 0.5;
    };

    auto niceEfficiency = [](const dDecay_t& event) { return invariantMass({event.kParams, event.pi1Params}) / 2; };

    // A less-nice efficiency that doesn't get recovered as nicely
    auto awkwardEfficiency = [](const dDecay_t& event) {
        std::vector<double> params = parametrisation(event);
        double              e{1};
        for (int i = 0; i < 3; ++i) {
            e *= params[i] / 2;
        }
        return 5 * e;
    };

    // An efficiency on the total pT of the k
    auto pTEfficiency = [](const dDecay_t& event) { return std::sqrt(pT(event.kParams)); };

    // Cast them to void so we can get away with only using one
    (void)simpleEfficiency;
    (void)niceEfficiency;
    (void)awkwardEfficiency;
    (void)pTEfficiency;

    // Create a thing for doing our efficiency correction
    size_t                    numBins    = 100;
    std::pair<double, double> axisLimits = std::make_pair(0., 2.);
    std::vector<double>       binLimits(numBins);
    for (size_t i = 0; i <= numBins; ++i) {
        binLimits[i] = axisLimits.first + (axisLimits.second - axisLimits.first) * i / (numBins);
    }

    PhspBins          Bins(5, binLimits);
    ChowLiuEfficiency EfficiencyCorrection(Bins);

    // Create a pseudo-Dalitz plot for the generated events
    // dalitzPlot(RootData.events, axisLimits, "All events", "ampgen.png");

    // Run the rejection thing on the AmpGen data
    std::vector<dDecay_t> detectedEvents = RootData.events;
    std::cout << "Run rejection thing on the AmpGen data..." << std::flush;
    applyEfficiency(&generator, pTEfficiency, detectedEvents);
    std::cout << "done" << std::endl;

    // Add the events of both type
    std::cout << "Add truth events..." << std::flush;
    for (auto truthEvent : RootData.events) {
        EfficiencyCorrection.addGeneratedEvent(parametrisation(truthEvent));
    }
    std::cout << "done" << std::endl;
    std::cout << "add detected events..." << std::flush;
    for (auto detectedEvent : detectedEvents) {
        EfficiencyCorrection.addMCEvent(parametrisation(detectedEvent));
    }
    std::cout << "done" << std::endl;

    // Perform efficiency parametrisation
    std::cout << "Perform efficiency parametrisation..." << std::flush;
    EfficiencyCorrection.efficiencyParametrisation();
    std::cout << "done" << std::endl;

    // Plot Daliz plot and histograms again
    // dalitzPlot(detectedEvents, axisLimits, "Events with efficiency", "eff.png");
    // histograms(RootData.events, detectedEvents, axisLimits, "Events with efficiency", "effAmpgen.png");

    correctionPlot(0, binLimits, RootData.events, detectedEvents, EfficiencyCorrection, "m12corr.png");
    correctionPlot(1, binLimits, RootData.events, detectedEvents, EfficiencyCorrection, "m23corr.png");
    correctionPlot(2, binLimits, RootData.events, detectedEvents, EfficiencyCorrection, "m34corr.png");
    correctionPlot(3, binLimits, RootData.events, detectedEvents, EfficiencyCorrection, "m123corr.png");
    correctionPlot(4, binLimits, RootData.events, detectedEvents, EfficiencyCorrection, "m234corr.png");

    return 0;
}
