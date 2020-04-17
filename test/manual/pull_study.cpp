#include <iostream>

#include "DecaySimulator.h"
#include "MinuitPolynomialFitter.h"
#include "RatioCalculator.h"

#include "FitterUtils.h"
#include "PullStudyHelpers.h"
#include "util.h"

#include <boost/progress.hpp>

void pull_study(const size_t meanNumCfEvents, const size_t numExperiments, bool integrate)
{
    // Choose parameters to use when simulating
    DecayParams_t phaseSpaceParams = {
        .x     = 0.0037,
        .y     = 0.0066,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = 2439.0,
    };
    double              maxTime           = 0.002;
    std::vector<double> expectedFitParams = util::expectedParams(phaseSpaceParams);
    size_t              numBins           = 50;

    // Create RNGs for numbers of decays
    double       meanNumDcsDecays = PullStudyHelpers::numDCSDecays(meanNumCfEvents, phaseSpaceParams, maxTime);
    std::mt19937 gen;
    std::poisson_distribution<size_t> cfDist(meanNumCfEvents);
    std::poisson_distribution<size_t> dcsDist(meanNumDcsDecays);

    // Find exponentially-spaced time bin limits to use
    std::vector<double> binLimits = PullStudyHelpers::exponentialBinLimits(maxTime, phaseSpaceParams.width, numBins);

    // Create a decay simulator
    SimulatedDecays MyDecays(maxTime, phaseSpaceParams);

    // Initialise vectors of fit parameter pulls and chi squared
    std::vector<double> aPull(numExperiments, -1);
    std::vector<double> bPull(numExperiments, -1);
    std::vector<double> cPull(numExperiments, -1);
    std::vector<double> chiSquaredVals(numExperiments, -1);

    boost::progress_display showProgress(numExperiments);

    for (size_t i = 0; i < numExperiments; ++i) {

        // Find how many decays to simulate
        size_t numCfEvents  = cfDist(gen);
        size_t numDcsEvents = dcsDist(gen);

        // Simulate them
        MyDecays.findCfDecayTimes(numCfEvents);
        MyDecays.findDcsDecayTimes(numDcsEvents);

        // Time binning
        std::vector<size_t> cfCounts  = util::binVector(MyDecays.RSDecayTimes, binLimits);
        std::vector<size_t> dcsCounts = util::binVector(MyDecays.WSDecayTimes, binLimits);

        // Find the ratio
        RatioCalculator MyRatios(cfCounts, dcsCounts, binLimits);
        MyRatios.calculateRatios();

        // Fit the ratio
        FitData_t              MyFitData(MyRatios.binCentres, MyRatios.binWidths, MyRatios.ratio, MyRatios.error);
        MinuitPolynomialFitter MyFitter(MyFitData, binLimits, phaseSpaceParams.width, integrate);
        MyFitter.setPolynomialParams(expectedFitParams, std::vector<double>(3, 1));
        MyFitter.fit();

        // Store parameter and chi squared
        aPull[i] = (MyFitter.fitParams.fitParams[0] - expectedFitParams[0]) / MyFitter.fitParams.fitParamErrors[0];
        bPull[i] = (MyFitter.fitParams.fitParams[1] - expectedFitParams[1]) / MyFitter.fitParams.fitParamErrors[1];
        cPull[i] = (MyFitter.fitParams.fitParams[2] - expectedFitParams[2]) / MyFitter.fitParams.fitParamErrors[2];
        chiSquaredVals[i] = *MyFitter.statistic;

        ++showProgress;
    }

    // Output mean and std dev of pulls
    std::vector<std::pair<double, double>> stats{};
    stats.push_back(PullStudyHelpers::meanAndStdDev(aPull));
    stats.push_back(PullStudyHelpers::meanAndStdDev(bPull));
    stats.push_back(PullStudyHelpers::meanAndStdDev(cPull));

    for (auto pair = stats.begin(); pair != stats.end(); ++pair) {
        std::cout << "Pull:\t" << pair->first << "+-" << pair->second << std::endl;
    }
}

int main(int argc, char* argv[])
{
    // If -i or --integrate specified, integrate over the bins
    bool        integrate{false};
    std::string iStr         = "-i";
    std::string integrateStr = "--integrate";
    for (int arg = 1; arg < argc; ++arg) {
        if (argv[arg] == iStr || argv[arg] == integrateStr) {
            integrate = true;
        }
    }

    std::string withIntegrating = integrate ? "with" : "without";
    std::cout << "Running pull study " << withIntegrating << " integrating over bins in fit." << std::endl;

    pull_study(10e5, 100, integrate);
}
