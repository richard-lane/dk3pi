/*
 * Perform a pull study for our fitter
 */

#include <algorithm>
#include <boost/progress.hpp>
#include <cassert>
#include <cmath>
#include <random>
#include <vector>

#include "TH1D.h"

#include "../lib/D2K3PiError.h"
#include "../lib/DecaySimulator.h"
#include "../lib/Fitter.h"
#include "../lib/PhaseSpaceBinning.h"
#include "../lib/RatioCalculator.h"
#include "PullStudyHelpers.h"
#include "util.h"

/*
 * Plot a histogram from a vector
 */
void plotHist(const std::vector<double>& vector, const size_t numBins, const std::string& name)
{
    assert((numBins != 0));
    double min = *(std::min_element(vector.begin(), vector.end()));
    double max = *(std::max_element(vector.begin(), vector.end()));

    double              binWidth = (max - min) / numBins;
    std::vector<double> binLimits(numBins + 1, -1);

    binLimits[0]       = min * 0.99;
    binLimits[numBins] = max * 1.01;
    for (size_t i = 1; i < numBins + 1; i++) {
        binLimits[i] = binWidth * i + min;
    }

    TH1D* hist = new TH1D(name.c_str(), name.c_str(), numBins, binLimits.data());
    hist->FillN(vector.size(), vector.data(), 0);

    util::saveObjectToFile(hist, name + ".pdf");
    delete hist;
}

/*
 * Perform a pull study with a specified number of experiments and events
 */
void pull_study(size_t nExperiments = 1000, size_t nEvents = 10000, size_t numPointsPerBin = 50)
{
    // Choose which parameters to use when simulating
    // These numbers are vaguely realistic but also entirely made up
    DecayParams_t phaseSpaceParams = {
        .x     = 0.0037,
        .y     = 0.0066,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = 2439.0,
    };

    // Restrict our times and rates to consider for our Monte Carlo simulation
    double maxTime = 0.002;

    // Calculate what we expect our fit parameters to be
    std::vector<double> expected_fit_params = util::expectedParams(phaseSpaceParams);
    double              expected_a          = expected_fit_params[0];
    double              expected_b          = expected_fit_params[1];
    double              expected_c          = expected_fit_params[2];

    // Create vectors of a, b and c
    std::vector<double> a_fit(nExperiments, -1);
    std::vector<double> b_fit(nExperiments, -1);
    std::vector<double> c_fit(nExperiments, -1);

    // Create vector to store our chi squared values
    std::vector<double> chiSqVector(nExperiments, -1);

    // Random number generators for finding how many of each type of decay to simulate
    double       meanNumDcsDecays = PullStudyHelpers::numDCSDecays(nEvents, phaseSpaceParams, maxTime);
    std::mt19937 generator;
    std::poisson_distribution<size_t> cfDistribution(nEvents);
    std::poisson_distribution<size_t> dcsDistribution(meanNumDcsDecays);

    // a quick and dirty progress bar
    // Note: this comes from a deprecated header, but no alterative exists yet
    boost::progress_display show_progress(nExperiments);

    // Store bin limits
    // This will be calculated on the first iteration
    std::vector<double> binLimits{};

    // Create decay simulator object
    SimulatedDecays MyDecays = SimulatedDecays(maxTime, phaseSpaceParams);

    // Initial guesses at fit parameters
    std::vector<double> parameterGuess{expected_a, expected_b, expected_c};
    std::vector<double> errorGuess{1, 1, 1};

    for (size_t i = 0; i < nExperiments; ++i) {
        // Find how many events of each type we need to generate
        // The number of events will have a well-defined mean, but will be drawn from Poisson distributions.
        size_t numCfEvents  = cfDistribution.operator()(generator);
        size_t numDcsEvents = dcsDistribution.operator()(generator);

        // Simulate DCS and CF decays with our parameters
        MyDecays.findDcsDecayTimes(numDcsEvents);
        MyDecays.findCfDecayTimes(numCfEvents);

        if (binLimits.empty()) {
            std::vector<double> dcsTimes{MyDecays.WSDecayTimes};
            std::sort(dcsTimes.begin(), dcsTimes.end());
            binLimits = util::findBinLimits(dcsTimes, numPointsPerBin, 0, 1.05 * maxTime);
        }

        // Plot histograms of event counts for both event types in each time bin
        // These will get saved as WSHist.pdf and RSHist.pdf
        // MyDecays.plotRates(timeBinLimits);

        RatioCalculator MyRatios = RatioCalculator(MyDecays.RSDecayTimes, MyDecays.WSDecayTimes, binLimits);
        MyRatios.calculateRatios();

        // Save the number of points stored in each bin to a text file
        // MyRatios.findNumPointsPerBin("numpoints.txt");

        // Fit our decays
        FitData_t MyFitData = FitData(MyRatios.binCentres, MyRatios.binWidths, MyRatios.ratio, MyRatios.error);
        Fitter    MyFitter  = Fitter(MyFitData);
        // MyFitter.fitUsingRootCustomFcn(0, maxTime * 1.2, "Q");
        MyFitter.fitUsingMinuit(parameterGuess, errorGuess, ChiSquared);

        // Save our fit plot to file
        if (!i) {
            const util::LegendParams_t legend = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = ""};
            std::string                path   = "fitplot_" + std::to_string(i) + ".pdf";
            MyFitter.saveFitPlot("Example fit plot", path, &legend);
        }

        // Store the parameters a, b and c
        // We care about their distance from the expected value, adjusted by their error
        a_fit[i] = (MyFitter.fitParams.fitParams[0] - expected_a) / MyFitter.fitParams.fitParamErrors[0];
        b_fit[i] = (MyFitter.fitParams.fitParams[1] - expected_b) / MyFitter.fitParams.fitParamErrors[1];
        c_fit[i] = (MyFitter.fitParams.fitParams[2] - expected_c) / MyFitter.fitParams.fitParamErrors[2];

        // Store the chi squared value in our vector
        chiSqVector[i] = *(MyFitter.statistic);

        ++show_progress;
    }

    // For each of a, b and c; plot the distance from the expected value
    PullStudyHelpers::plot_parameter_distribution("a", a_fit, nExperiments);
    PullStudyHelpers::plot_parameter_distribution("b", b_fit, nExperiments);
    PullStudyHelpers::plot_parameter_distribution("c", c_fit, nExperiments);

    // Plot the distribution of chi squared values
    plotHist(chiSqVector, 50, "MinuitChiSq");
}

// Hide this program's main() function from ROOT's Cling interpreter
#ifndef __CINT__
int main()
{
    pull_study(100, 800000, 50);

    return 0;
}
#endif // __CINT__
