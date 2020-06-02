/*
 * Run a simulation using a SimulatedDecays object and draw a graph of the two curves.
 *
 * i wrote this during lhcb week so it is hangover-quality code
 */
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <utility>

#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>

#include "../pull_study/DecaySimulator.h"
#include "../pull_study/PullStudyHelpers.h"
#include "RatioCalculator.h"
#include "physics.h"
#include "util.h"

/*
 * Find chi squared between a dataset with errors and a TF1
 */
double chiSq(const std::vector<double> &times,
             const std::vector<double> &data,
             const std::vector<double> &errors,
             const TF1 *                fcn)
{

    double chi2 = 0.0;
    for (size_t i = 0; i < times.size(); ++i) {
        chi2 += std::pow((fcn->Eval(times[i]) - data[i]) / errors[i], 2);
    }
    return chi2;
}

/*
 * Using the forms of the functions given in S Harnew's paper, generate a histogram of expected of DCS and CF events
 * The, use the decay simulator to generate the right numbers of DCS and CF events (we want to have the same numbers of
 * analytical + accept-reject events for both decay types)
 *
 * Plot two histograms with acc-rej and analytical histograms overlaid, WS.pdf and RS.pdf
 *
 * This function Sucks but its a script that I only really intend to run once, so thats probably ok
 */
void simulateDecays()
{
    // Parameter values
    DecayParams_t MyParams{
        .x     = 0.004,
        .y     = 0.007,
        .r     = 0.055,
        .z_im  = -0.3,
        .z_re  = 0.7,
        .width = 2500.0,
    };
    double maxTime = 0.002;

    // Define time bin limits
    size_t              numTimeBins = 50;
    std::vector<double> timeBinLimits{};
    for (size_t i = 0; i < numTimeBins + 1; ++i) {
        timeBinLimits.push_back(i * (maxTime / numTimeBins));
    }

    // Choose about how many decays we want
    size_t numDecays = 10000000;

    // Find how many events we expect in each bin
    double efficiencyTimescale = 1 / MyParams.width;
    auto   rsRate              = [&](const double x) { return Phys::cfRate(x, MyParams, efficiencyTimescale); };
    auto   wsRate              = [&](const double x) { return Phys::dcsRate(x, MyParams, efficiencyTimescale); };
    double cfIntegral          = util::gaussLegendreQuad(rsRate, 0, maxTime);
    double dcsIntegral         = util::gaussLegendreQuad(wsRate, 0, maxTime);

    std::vector<double> expectedCfBinPopulation(numTimeBins, -1);
    std::vector<double> expectedDcsBinPopulation(numTimeBins, -1);
    for (size_t i = 0; i < numTimeBins; ++i) {
        expectedCfBinPopulation[i] =
            numDecays * util::gaussLegendreQuad(rsRate, timeBinLimits[i], timeBinLimits[i + 1]) / cfIntegral;
        expectedDcsBinPopulation[i] =
            numDecays * util::gaussLegendreQuad(wsRate, timeBinLimits[i], timeBinLimits[i + 1]) / dcsIntegral;
    }

    // Generator and PDF for random numbers
    std::random_device                     rd;
    std::shared_ptr<std::mt19937>          _gen = std::make_shared<std::mt19937>(rd());
    std::uniform_real_distribution<double> uniform;

    auto gen = [&](void) {
        double x = uniform(*_gen);
        double z = 1 - std::exp(-1 * MyParams.width * maxTime);
        return (-1 / MyParams.width) * std::log(1 - z * x);
    };
    auto genPDF = [&](double x) {
        return std::exp(-MyParams.width * x) * MyParams.width / (1 - std::exp(-MyParams.width * maxTime));
    };

    SimulatedDecays MyDecays = SimulatedDecays(gen, genPDF, rsRate, wsRate, std::make_pair(0., maxTime), _gen);
    MyDecays.findCfDecayTimes(numDecays);
    MyDecays.findDcsDecayTimes(numDecays);

    // Test the PDF and generator give sensible stuff
    MyDecays.test(1e6, timeBinLimits);

    // Plot both the expected and Monte-Carlo generated decay rates on the same axes
    TH1D *generatedRSHist = new TH1D("Test accept-reject, RS",
                                     "RS expected and Monte-Carlo event numbers;time/ns",
                                     numTimeBins - 1,
                                     timeBinLimits.data());
    TH1D *generatedWSHist = new TH1D("Test accept-reject, WS",
                                     "WS expected and Monte-Carlo event numbers;time/ns",
                                     numTimeBins - 1,
                                     timeBinLimits.data());
    TH1D *expectedRSHist  = new TH1D("expected, RS", "", numTimeBins - 1, timeBinLimits.data());
    TH1D *expectedWSHist  = new TH1D("expected, WS", "", numTimeBins - 1, timeBinLimits.data());

    // Histogram of our acc-rej data
    generatedRSHist->FillN(numDecays, MyDecays.RSDecayTimes.data(), nullptr);
    generatedWSHist->FillN(numDecays, MyDecays.WSDecayTimes.data(), nullptr);

    // Histogram of our expected data
    // Bin numbering starts at 1 for some reason
    for (size_t i = 1; i <= numTimeBins; ++i) {
        expectedRSHist->SetBinContent(i, expectedCfBinPopulation[i - 1]);
        expectedWSHist->SetBinContent(i, expectedDcsBinPopulation[i - 1]);
    }

    TCanvas *rsCanvas = new TCanvas();
    rsCanvas->cd();
    generatedRSHist->Draw("");
    expectedRSHist->Draw("CSAME");
    rsCanvas->SaveAs("RS.pdf");
    delete rsCanvas;

    TCanvas *wsCanvas = new TCanvas();
    wsCanvas->cd();
    generatedWSHist->Draw("");
    expectedWSHist->Draw("CSAME");
    wsCanvas->SaveAs("WS.pdf");
    delete wsCanvas;

    // Take the ratio of our datasets, plot against the expected ratio
    std::vector<double> expectedParams    = util::expectedParams(MyParams);
    TF1 *               expectedRatioFunc = new TF1("expectedPoly", "[0] + [1]*x + [2]*x*x", 0, maxTime);
    expectedRatioFunc->SetParameter(0, expectedParams[0]);
    expectedRatioFunc->SetParameter(1, expectedParams[1]);
    expectedRatioFunc->SetParameter(2, expectedParams[2]);

    // Recalculate the number of DCS events using the right scaling
    MyDecays.findDcsDecayTimes(PullStudyHelpers::numDCSDecays(numDecays, MyParams, maxTime, efficiencyTimescale));
    std::vector<size_t> cfCounts  = util::binVector(MyDecays.RSDecayTimes, timeBinLimits);
    std::vector<size_t> dcsCounts = util::binVector(MyDecays.WSDecayTimes, timeBinLimits);

    RatioCalculator Calculator = RatioCalculator(cfCounts, dcsCounts, timeBinLimits);
    Calculator.calculateRatios();

    TGraphErrors *ratioGraph = new TGraphErrors(
        numTimeBins, Calculator.binCentres.data(), Calculator.ratio.data(), nullptr, Calculator.error.data());
    ratioGraph->SetTitle("Expected vs Actual DCS/CF ratio;time/ns;DCS/CF counts per bin");

    // Find chi sq
    double chiSquare = chiSq(Calculator.binCentres, Calculator.ratio, Calculator.error, expectedRatioFunc);
    std::cout << "numbins :" << numTimeBins << std::endl;
    std::cout << "chiSquare: " << chiSquare << std::endl;

    TCanvas *ratioCanvas = new TCanvas();
    ratioCanvas->cd();
    ratioGraph->Draw("AP");
    expectedRatioFunc->Draw("CSAME");
    wsCanvas->SaveAs("ratio.pdf");
    delete ratioCanvas;

    delete generatedRSHist;
    delete generatedWSHist;
    delete expectedRSHist;
    delete expectedWSHist;
    delete ratioGraph;
    delete expectedRatioFunc;
}

int main()
{
    simulateDecays();

    return 0;
}
