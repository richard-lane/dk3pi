/*
 * Run a simulation using a SimulatedDecays object and draw a graph of the two curves.
 *
 * i wrote this during lhcb week so it is hangover-quality code
 */
#include <iostream>
#include <numeric>
#include <random>
#include <utility>

#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>

#include "DecaySimulator.h"
#include "MinuitFitter.h"
#include "PullStudyHelpers.h"
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
 * Set the decay parameters to sensible values.
 */
void setParams(DecayParams_t &DecayParams)
{
    // Global Params
    DecayParams.x     = 0.004;
    DecayParams.y     = 0.007;
    DecayParams.width = 2500.0; // Width is in nanoseconds

    // Per-bin params
    DecayParams.r    = 1.14;
    DecayParams.z_re = 0.7;
    DecayParams.z_im = -0.3;
}

/*
 * Takes in an approximate number of points, some bin limits and a function pointer (DecayParams_t, double -> double);
 * uses them to generate a distribution of points that describes the function
 *
 * Uses 1st order trapezium rule (i.e. approximates the area of each bin with a trapezium)
 *
 */
std::vector<size_t> expectedFunction(const size_t               approxNumPoints,
                                     const DecayParams_t &      MyParams,
                                     const std::vector<double> &binLimits,
                                     double (*func)(const double, const DecayParams_t &))
{
    // Evaluate func(time) * bigNumber to get a large number in each bin
    // Otherwise wen we cast our double func() value to a size_t it will just end up being 0
    size_t              bigNumber = 1000000000;
    size_t              numBins   = binLimits.size() - 1;
    std::vector<size_t> values(numBins, (size_t)NAN);

    for (size_t i = 0; i < numBins; ++i) {
        double binWidth  = binLimits[i + 1] - binLimits[i];
        double lowerFunc = func(binLimits[i], MyParams);
        double upperFunc = func(binLimits[i + 1], MyParams);
        values[i]        = bigNumber * 0.5 * binWidth * (lowerFunc + upperFunc);
    }

    // Want to rescale our expected values such that we have numDecays decays total
    size_t numEvents = std::accumulate(values.begin(), values.end(), (size_t)0);
    std::transform(
        values.begin(),
        values.end(),
        values.begin(),
        std::bind(std::multiplies<double>(), std::placeholders::_1, (double)approxNumPoints / (double)numEvents));

    return values;
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
    DecayParams_t MyParams;
    setParams(MyParams);
    double maxTime = 0.002;

    // Define time bin limits
    size_t              numTimeBins = 50;
    std::vector<double> timeBinLimits{};
    for (size_t i = 0; i < numTimeBins + 1; ++i) {
        timeBinLimits.push_back((1 / MyParams.width) + i * (maxTime - 1 / MyParams.width) / numTimeBins);
    }

    // Choose about how many decays we want
    // This won't be exact- we convert a vector of doubles to a vector of size_t, so we end don't get the exact number
    // we want
    size_t approxNumDecays = 100000;

    // Take the centre of each time bin, calculate the rate in each and multiply it by a large number to get many events
    // in each bin
    std::vector<size_t> expectedCfBinPopulation =
        expectedFunction(approxNumDecays, MyParams, timeBinLimits, Phys::rightSignDecayRate);
    std::vector<size_t> expectedDcsBinPopulation =
        expectedFunction(approxNumDecays, MyParams, timeBinLimits, Phys::wrongSignDecayRate);

    // After rescaling, count how many of each event type we have (should be about approxNumDecays)
    size_t numCf  = std::accumulate(expectedCfBinPopulation.begin(), expectedCfBinPopulation.end(), (size_t)0);
    size_t numDcs = std::accumulate(expectedDcsBinPopulation.begin(), expectedDcsBinPopulation.end(), (size_t)0);

    // Create our Decay simulator object and generate CF and DCS times using accept-reject
    SimulatedDecays MyDecays = SimulatedDecays(maxTime, MyParams);
    MyDecays.findCfDecayTimes(numCf);
    MyDecays.findDcsDecayTimes(numDcs);

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
    generatedRSHist->FillN(numCf, MyDecays.RSDecayTimes.data(), nullptr);
    generatedWSHist->FillN(numDcs, MyDecays.WSDecayTimes.data(), nullptr);

    // Histogram of our expected data
    // Bin numbering starts at 1 for some reason
    for (size_t i = 1; i <= numTimeBins; ++i) {
        expectedRSHist->SetBinContent(i, expectedCfBinPopulation[i - 1]);
        expectedWSHist->SetBinContent(i, expectedDcsBinPopulation[i - 1]);
    }

    TCanvas *rsCanvas = new TCanvas();
    rsCanvas->cd();
    generatedRSHist->Draw("");
    generatedRSHist->SetMinimum(0);
    generatedRSHist->SetMaximum(numCf / (8 * numTimeBins / 50));
    expectedRSHist->Draw("CSAME");
    rsCanvas->SaveAs("RS.pdf");
    delete rsCanvas;

    TCanvas *wsCanvas = new TCanvas();
    wsCanvas->cd();
    generatedWSHist->Draw("");
    generatedWSHist->SetMinimum(0);
    generatedWSHist->SetMaximum(
        numCf / (8 * numTimeBins / 50)); // Use numCf to set axis limits so they are the same for both histograms
    expectedWSHist->Draw("CSAME");
    wsCanvas->SaveAs("WS.pdf");
    delete wsCanvas;

    // Take the ratio of our datasets, plot against the expected ratio
    std::vector<double> expectedParams    = util::expectedParams(MyParams);
    TF1 *               expectedRatioFunc = new TF1("expectedPoly", "[0] + [1]*x + [2]*x*x", 0, maxTime);
    expectedRatioFunc->SetParameter(0, expectedParams[0]);
    expectedRatioFunc->SetParameter(1, expectedParams[1]);
    expectedRatioFunc->SetParameter(2, expectedParams[2]);

    numDcs = (size_t)PullStudyHelpers::numDCSDecays(numCf, MyParams, maxTime);
    MyDecays.findDcsDecayTimes(numDcs);

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
