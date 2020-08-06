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

#include "ReadRoot.h"
#include "efficiencyUtil.h"
#include "toyStudy.h"
#include "util.h"

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

void histograms(const std::vector<dDecay_t>      events,
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

    size_t                numBins = 100;
    std::unique_ptr<TH1D> m12Hist(
        new TH1D("m12", (title + " m12").c_str(), numBins, axisLimits.first, axisLimits.second));
    std::unique_ptr<TH1D> m13Hist(
        new TH1D("m23", (title + " m13").c_str(), numBins, axisLimits.first, axisLimits.second));

    m12Hist->FillN(numEvents, m12.data(), nullptr);
    m13Hist->FillN(numEvents, m13.data(), nullptr);

    std::unique_ptr<TCanvas> c(new TCanvas());
    c->Divide(1, 2);
    c->cd(1);
    m12Hist->Draw();
    c->cd(2);
    m13Hist->Draw();

    c->SaveAs(path.c_str());
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
    ReadRoot                 RootData(tFile.get(), ampgenTreeName, ampgenBranchNames, ampgenMomentumPostfixes);

    std::pair<double, double> axisLimits = std::make_pair(0.5, 1.8);
    dalitzPlot(RootData.events, axisLimits, "All events", "ampgen.png");
    histograms(RootData.events, axisLimits, "All events", "ampgenHists.png");

    // Create a lambda fcn that will reject all events with m12 below some value that i can read off the Dalitz plot
    // when it comes to it
    // Let's say that events with m12 > 1 aren't detected, but all others are
    auto detectionChance = [&](const dDecay_t& event) {
        if (invariantMass(std::vector<kinematicParams_t>{event.kParams, event.pi1Params}) < 1) {
            return 1.;
        }
        return 0.5;
    };

    // Run the rejection thing on the AmpGen data
    std::random_device rd;
    std::mt19937       generator(rd());
    applyEfficiency(&generator, detectionChance, RootData.events);

    // Plot Daliz plot and histograms again
    dalitzPlot(RootData.events, axisLimits, "Events with efficiency", "eff.png");
    histograms(RootData.events, axisLimits, "Events with efficiency", "effAmpgen.png");

    return 0;
}
