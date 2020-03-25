/*
 * Perform a pull study for our fitter
 */

#include <algorithm>
#include <boost/progress.hpp>
#include <random>
#include <vector>

#include "TH1D.h"

#include "../lib/DecaySimulator.h"
#include "../lib/Fitter.h"
#include "../lib/MCGenerator.h"
#include "../lib/PhaseSpaceBinning.h"
#include "../lib/RatioCalculator.h"
#include "PullStudyHelpers.h"
#include "util.h"

/*
 * Perform a pull study with a specified number of experiments and events
 */
void pull_study(size_t nExperiments = 1000, size_t nEvents = 10000)
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
    double                    maxTime      = 0.002;
    std::pair<double, double> allowedTimes = std::make_pair(0, maxTime);
    std::pair<double, double> allowedRates = std::make_pair(0, 1.3);

    // Create sensible time bins
    // size_t              numTimeBins = 500;
    // std::vector<double> timeBinLimits{};
    // for (size_t i = 0; i < numTimeBins; ++i) {
    //    timeBinLimits.push_back(i * maxTime * 1.1 / numTimeBins);
    //}

    // Calculate what we expect our fit parameters to be
    std::vector<double> expected_fit_params = PullStudyHelpers::expectedParams(phaseSpaceParams);
    double              expected_a          = expected_fit_params[0];
    double              expected_b          = expected_fit_params[1];
    double              expected_c          = expected_fit_params[2];

    // Create vectors of a, b and c
    std::vector<double> a_fit(nExperiments, -1);
    std::vector<double> b_fit(nExperiments, -1);
    std::vector<double> c_fit(nExperiments, -1);

    // Random number generators for finding how many of each type of decay to simulate
    size_t       meanNumDcsDecays = PullStudyHelpers::numDCSDecays(nEvents, phaseSpaceParams, maxTime);
    std::mt19937 generator;
    std::poisson_distribution<size_t> cfDistribution(nEvents);
    std::poisson_distribution<size_t> dcsDistribution(meanNumDcsDecays);

    // a quick and dirty progress bar
    // Note: this comes from a deprecated header, but no alterative exists yet
    boost::progress_display show_progress(nExperiments);

    for (size_t i = 0; i < nExperiments; ++i) {
        // Find how many events of each type we need to generate
        // The number of events will have a well-defined mean, but will be drawn from Poisson distributions.
        size_t numCfEvents  = cfDistribution.operator()(generator);
        size_t numDcsEvents = dcsDistribution.operator()(generator);

        // Simulate DCS and CF decays with our parameters
        SimulatedDecays MyDecays = SimulatedDecays(allowedTimes, allowedRates, phaseSpaceParams);
        MyDecays.findDcsDecayTimes(numDcsEvents);
        MyDecays.findCfDecayTimes(numCfEvents);

        // Find bin limits such that we have at least 50 DCS points per bin
        std::vector<double> dcsTimes{MyDecays.WSDecayTimes};
        std::sort(dcsTimes.begin(), dcsTimes.end());
        std::vector<double> timeBinLimits = util::findBinLimits(dcsTimes, 50, 0, 1.05 * maxTime);

        // Plot histograms of event counts for both event types in each time bin
        // These will get saved as WSHist.pdf and RSHist.pdf
        // MyDecays.plotRates(timeBinLimits);

        RatioCalculator MyRatios = RatioCalculator(MyDecays.RSDecayTimes, MyDecays.WSDecayTimes, timeBinLimits);
        MyRatios.calculateRatios();

        // Save the number of points stored in each bin to a text file
        MyRatios.findNumPointsPerBin("numpoints.txt");

        // Fit our decays
        FitData_t MyFitData = FitData(MyRatios.binCentres, MyRatios.binWidths, MyRatios.ratio, MyRatios.error);
        Fitter    MyFitter  = Fitter(MyFitData);
        // MyFitter.fitUsingRootCustomFcn(0, maxTime * 1.2, "Q");

        // Initially guess the parameters are their known values...
        // FitMethod is either ChiSquared or MaxLikelihood
        std::vector<double> parameterGuess{expected_a, expected_b, expected_c};
        std::vector<double> errorGuess{1, 1, 1};
        MyFitter.fitUsingMinuit(parameterGuess, errorGuess, ChiSquared);

        // Save our fit plot to file
        // std::string path = "fitplot_" + std::to_string(i) + ".pdf";
        // MyFitter.saveFitPlot("foo", path);

        // Store the parameters a, b and c
        // We care about their distance from the expected value, adjusted by their error
        a_fit[i] = (MyFitter.fitParams.fitParams[0] - expected_a) / MyFitter.fitParams.fitParamErrors[0];
        b_fit[i] = (MyFitter.fitParams.fitParams[1] - expected_b) / MyFitter.fitParams.fitParamErrors[1];
        c_fit[i] = (MyFitter.fitParams.fitParams[2] - expected_c) / MyFitter.fitParams.fitParamErrors[2];

        ++show_progress;
    }

    // For each of a, b and c; plot the distance from the expected value
    PullStudyHelpers::plot_parameter_distribution("a", a_fit, nExperiments);
    PullStudyHelpers::plot_parameter_distribution("b", b_fit, nExperiments);
    PullStudyHelpers::plot_parameter_distribution("c", c_fit, nExperiments);
}

// Hide this program's main() function from ROOT's Cling interpreter
#ifndef __CINT__
int main()
{
    pull_study(250, 100000);

    return 0;
}
#endif // __CINT__
