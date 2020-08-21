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
#include "util.h"

struct EventDetectionProbNotNormalised : public std::exception {
    EventDetectionProbNotNormalised(const double prob)
        : _msg("Probability " + std::to_string(prob) + " not between 0 and 1")
    {
        ;
    }

    const char* what() const throw() { return _msg.c_str(); }

  private:
    const std::string _msg;
};

/*
 * Provide a random number generator, a vector of D decay events and a function that provides an event's detection
 * probability
 *
 * The probability of detecting an event should be a number between 0 and 1
 *
 * Removes the undetected events from the vector in-place.
 */
void applyEfficiency(std::mt19937* const                    generator,
                     const std::function<double(dDecay_t)>& eventDetectionProb,
                     std::vector<dDecay_t>&                 events)
{
    std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);

    // Lambda that we'll use to decide whether to remove our event
    auto removeEvent = [&](const dDecay_t& event) {
        double detectionProb = eventDetectionProb(event);
        if (detectionProb < 0 || detectionProb > 1) {
            throw EventDetectionProbNotNormalised(detectionProb);
        }
        return uniformDistribution(*generator) > detectionProb;
    };

    // Move the items to remove to the end of the vector and erase them
    auto it = std::remove_if(events.begin(), events.end(), removeEvent);
    events.erase(it, events.end());
}

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

double simpleEfficiency(const dDecay_t& event)
{
    (void)event;
    return 0.5;
}

/*
 * A nice efficiency that gets recovered well
 */
double niceEfficiency(const dDecay_t& event)
{
    return invariantMass({event.kParams, event.pi1Params}) / 2;
}

/*
 * A less-nice efficiency that doesn't get recovered as nicely
 */
double awkwardEfficiency(const dDecay_t& event)
{
    std::vector<double> params = parametrisation(event);
    double              e{1};
    for (int i = 0; i < 3; ++i) {
        e *= params[i] / 2;
    }
    return 5 * e;
}

/*
 * An efficiency on the total pT of the k
 */
double pTEfficiency(const dDecay_t& event)
{
    return std::sqrt(pT(event.kParams));
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

    // Create a thing for doing our efficiency correction
    size_t                    numBins    = 100;
    std::pair<double, double> axisLimits = std::make_pair(0., 2.);
    std::vector<double>       binLimits(numBins);
    for (size_t i = 0; i <= numBins; ++i) {
        binLimits[i] = axisLimits.first + (axisLimits.second - axisLimits.first) * i / (numBins);
    }

    PhspBins          Bins(5, binLimits);
    ChowLiuEfficiency EfficiencyCorrection(Bins);

    // Run the rejection thing on the AmpGen data
    std::cout << "Run rejection on the AmpGen data..." << std::flush;
    std::vector<dDecay_t> detectedEvents = RootData.events;
    std::random_device    rd;
    std::mt19937          generator(rd());
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

    correctionPlot(0, binLimits, RootData.events, detectedEvents, EfficiencyCorrection, "m12corr.png");
    correctionPlot(1, binLimits, RootData.events, detectedEvents, EfficiencyCorrection, "m23corr.png");
    correctionPlot(2, binLimits, RootData.events, detectedEvents, EfficiencyCorrection, "m34corr.png");
    correctionPlot(3, binLimits, RootData.events, detectedEvents, EfficiencyCorrection, "m123corr.png");
    correctionPlot(4, binLimits, RootData.events, detectedEvents, EfficiencyCorrection, "m234corr.png");

    return 0;
}
