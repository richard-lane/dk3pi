/*
 * Perform a pull study for our fitter
 */

#include <iostream>
#include <vector>

#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TRandom.h"

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
 */
void plot_parameter_distribution(std::string title, std::vector<double> parameter, size_t nExperiments)
{
    TH1D *MyGraph = new TH1D(title.c_str(), title.c_str(), 200, -2, 2);
    MyGraph->FillN(nExperiments, parameter.data(), 0);
    util::saveToFile(MyGraph, (title + ".pdf").c_str());
    delete MyGraph;
}

/*
 * Perform a pull study with a specified number of experiments and events
 */
void pull_study(size_t nExperiments = 1000, size_t nEvents = 10000)
{
    // Choose which parameters to use when simulating
    DecayParams_t phaseSpaceParams = {
        .x     = 0.0037,
        .y     = 0.0066,
        .r     = 1.1319,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = 2439.0,
    };

    // Restrict our times and rates to consider for our Monte Carlo simulation
    double                    maxTime      = 0.002;
    std::pair<double, double> allowedTimes = std::make_pair(0, maxTime);
    std::pair<double, double> allowedRates = std::make_pair(0, 1.3);

    // Create sensible time bins
    size_t              numTimeBins = 200;
    std::vector<double> timeBinLimits{};
    for (size_t i = 0; i < numTimeBins; ++i) {
        timeBinLimits.push_back(i * i * i * maxTime * 1.1 / (numTimeBins * numTimeBins * numTimeBins));
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
        SimulatedDecays MyDecays = SimulatedDecays(allowedTimes, allowedRates, phaseSpaceParams, nEvents);
        size_t          nWS      = 0;
        size_t          nRS      = 0;
        // Keep generating points, checking if they are in either distribution until both the RS and WS vectors are
        // populated.
        // This case uses the same random generation for both RS and WS data- we probably want to do different random
        // numbers for each
        while (nWS < nEvents && nRS < nEvents) {
            double time  = MyDecays.getRandomX();
            double ratio = MyDecays.getRandomY();

            if (nWS < nEvents && MyDecays.isAccepted(time, ratio, false)) {
                MyDecays.WSDecayTimes[nWS] = time;
                nWS++;
            }

            if (nRS < nEvents && MyDecays.isAccepted(time, ratio, true)) {
                MyDecays.RSDecayTimes[nRS] = time;
                nRS++;
            }
        }

        // Calculate the ratio between our two decay time sets
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
    pull_study(5000, 10000);
}
#endif // __CINT__
