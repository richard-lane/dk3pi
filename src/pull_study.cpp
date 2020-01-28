/*
 * Perform a pull study for our fitter
 */

#include <iostream>
#include <vector>

#include "TFile.h"
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
 * Perform a pull study with a specified number of experiments and events
 */
void pull_study(size_t nExperiments = 100, size_t nEvents = 10000)
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
    std::cerr << nExperiments + nEvents;

    // Calculate what we expect our fit parameters to be

    for (size_t i = 0; i < nExperiments; ++i) {
        // Simulate DCS and CF decays with our parameters
        std::cerr << phaseSpaceParams.x;

        // Fit our decays

        // Store the parameters a, b and c
    }

    // For each of a, b and c; plot the distance from the expected value
}

// Hide this program's main() function from ROOT's Cling interpreter
#ifndef __CINT__
int main()
{
    pull_study();
}
#endif // __CINT__
