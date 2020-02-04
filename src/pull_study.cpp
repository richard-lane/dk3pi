/*
 * Perform a pull study for our fitter
 */

#include <vector>

#include "TH1D.h"

#include "util.cpp"

#include "../include/DecaySimulator.h"
#include "../include/Fitter.h"
#include "../include/MCGenerator.h"
#include "../include/PhaseSpaceBinning.h"
#include "../include/RatioCalculator.h"

#include "../src/DecaySimulator.cpp"
#include "../src/Fitter.cpp"
#include "../src/MCGenerator.cpp"
#include "../src/PhaseSpaceBinning.cpp"
#include "../src/RatioCalculator.cpp"

/*
 * Plot the distribution of a vector as a histogram of 200 bins between -2 and 2
 *
 * Plot is centred on expectedMean- this should be 0 but it isn't, so i'm doing this to make the plots sensible
 */
void plot_parameter_distribution(std::string         title,
                                 std::vector<double> parameter,
                                 size_t              nExperiments,
                                 double              expectedMean  = 0,
                                 double              expectedSigma = 1)
{
    // Define axis limits
    double xMin = expectedMean - 5 * expectedSigma;
    double xMax = expectedMean - 5 * expectedSigma;

    TH1D *MyGraph = new TH1D(title.c_str(), title.c_str(), 200, xMin, xMax);

    MyGraph->FillN(nExperiments, parameter.data(), 0);
    MyGraph->SetTitle((title + ";Normalised Pull;Count").c_str());

    util::saveToFile(MyGraph, (title + ".pdf").c_str());

    std::cout << title + " mean:\t\t" + MyGraph->GetMean() << std::endl;
    std::cout << title + " std dev:\t" + MyGraph->GetStdDev() << std::endl;
    delete MyGraph;
}

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

    size_t numCfEvents  = nEvents;
    size_t numDcsEvents = nEvents * phaseSpaceParams.r * phaseSpaceParams.r;

    // Restrict our times and rates to consider for our Monte Carlo simulation
    double                    maxTime      = 0.002;
    std::pair<double, double> allowedTimes = std::make_pair(0, maxTime);
    std::pair<double, double> allowedRates = std::make_pair(0, 1.3);

    // Create sensible time bins
    size_t              numTimeBins = 20;
    std::vector<double> timeBinLimits{};
    for (size_t i = 0; i < numTimeBins; ++i) {
        timeBinLimits.push_back(i * maxTime * 1.1 / numTimeBins);
    }

    // Calculate what we expect our fit parameters to be
    double expected_a = phaseSpaceParams.r * phaseSpaceParams.r;
    double expected_b = phaseSpaceParams.r *
                        (phaseSpaceParams.y * phaseSpaceParams.z_re + phaseSpaceParams.x * phaseSpaceParams.z_im) *
                        phaseSpaceParams.width;
    double expected_c = 0.25 * (std::pow(phaseSpaceParams.x, 2) + std::pow(phaseSpaceParams.y, 2)) *
                        std::pow(phaseSpaceParams.width, 2);

    // Create vectors of a, b and c
    std::vector<double> a_fit(nExperiments, -1);
    std::vector<double> b_fit(nExperiments, -1);
    std::vector<double> c_fit(nExperiments, -1);

    for (size_t i = 0; i < nExperiments; ++i) {
        // Simulate DCS and CF decays with our parameters
        SimulatedDecays MyDecays = SimulatedDecays(allowedTimes, allowedRates, phaseSpaceParams);

        MyDecays.findDcsDecayTimes(numDcsEvents);
        MyDecays.findCfDecayTimes(numCfEvents);

        RatioCalculator MyRatios = RatioCalculator(MyDecays.RSDecayTimes, MyDecays.WSDecayTimes, timeBinLimits);
        MyRatios.calculateRatios();

        // Fit our decays
        FitData_t MyFitData = FitData(MyRatios.binCentres, MyRatios.binWidths, MyRatios.ratio, MyRatios.error);
        Fitter    MyFitter  = Fitter(MyFitData);
        MyFitter.pol2fit("Q");

        // Store the parameters a, b and c
        // We care about their distance from the expected value, adjusted by their error
        a_fit[i] = (MyFitter.fitParams.fitParams[0] - expected_a) / MyFitter.fitParams.fitParamErrors[0];
        b_fit[i] = (MyFitter.fitParams.fitParams[1] - expected_b) / MyFitter.fitParams.fitParamErrors[1];
        c_fit[i] = (MyFitter.fitParams.fitParams[2] - expected_c) / MyFitter.fitParams.fitParamErrors[2];
    }

    // For each of a, b and c; plot the distance from the expected value
    plot_parameter_distribution("a", a_fit, nExperiments);
    plot_parameter_distribution("b", b_fit, nExperiments);
    plot_parameter_distribution("c", c_fit, nExperiments);
}

// Hide this program's main() function from ROOT's Cling interpreter
#ifndef __CINT__
int main()
{
    pull_study(2000, 100000);
}
#endif // __CINT__
