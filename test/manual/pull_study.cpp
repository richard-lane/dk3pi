#include <iostream>

#include "DecaySimulator.h"
#include "MinuitPolynomialFitter.h"
#include "PhysicalFitter.h"
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
    double maxTime = 10 / phaseSpaceParams.width;
    // std::vector<double> expectedFitParams = util::expectedParams(phaseSpaceParams);
    size_t numBins = 50;

    // Create RNGs for numbers of decays
    double       meanNumDcsDecays = PullStudyHelpers::numDCSDecays(meanNumCfEvents, phaseSpaceParams, maxTime);
    std::mt19937 gen;
    std::poisson_distribution<size_t> cfDist(meanNumCfEvents);
    std::poisson_distribution<size_t> dcsDist(meanNumDcsDecays);

    // Find exponentially-spaced time bin limits to use
    std::vector<double> binLimits =
        util::exponentialBinLimits(1 / phaseSpaceParams.width, maxTime, phaseSpaceParams.width, numBins);

    // Create a decay simulator
    SimulatedDecays MyDecays(maxTime, phaseSpaceParams);

    // Initialise vectors of fit parameter pulls and chi squared
    // std::vector<double> aPull(numExperiments, -1);
    // std::vector<double> bPull(numExperiments, -1);
    // std::vector<double> cPull(numExperiments, -1);
    std::vector<double> rPull(numExperiments, -1);
    std::vector<double> imZPull(numExperiments, -1);
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

        // Create integral options if we need them
        std::unique_ptr<IntegralOptions_t> integralOptions = nullptr;
        if (integrate) {
            integralOptions = std::make_unique<IntegralOptions_t>(phaseSpaceParams.width, binLimits, 1e-10, 10);
        }
        FitData_t MyFitData(MyRatios.binCentres, MyRatios.binWidths, MyRatios.ratio, MyRatios.error);

        // Fit data
        // MinuitPolynomialFitter MyFitter(MyFitData, integralOptions.get());
        // MyFitter.setPolynomialParams(expectedFitParams, std::vector<double>(3, 1));
        // MyFitter.fixParameters(std::vector<std::string>{"a"});
        // MyFitter.fit();

        PhysicalFitter MyFitter(MyFitData, *integralOptions);
        MyFitter.setPhysicalFitParams(std::vector<double>{phaseSpaceParams.x,
                                                          phaseSpaceParams.y,
                                                          phaseSpaceParams.r,
                                                          phaseSpaceParams.z_im,
                                                          phaseSpaceParams.z_re,
                                                          phaseSpaceParams.width},
                                      std::vector<double>(6, 1));
        MyFitter.fixParameters(std::vector<std::string>{"width", "z_re", "x", "y"});
        MyFitter.fit();

        // Store parameter and chi squared
        rPull[i]   = (MyFitter.fitParams.fitParams[2] - phaseSpaceParams.r) / MyFitter.fitParams.fitParamErrors[2];
        imZPull[i] = (MyFitter.fitParams.fitParams[3] - phaseSpaceParams.z_im) / MyFitter.fitParams.fitParamErrors[3];
        chiSquaredVals[i] = *MyFitter.statistic;

        ++showProgress;
    }

    // Output mean and std dev of pulls
    std::vector<std::pair<double, double>> stats{};
    // stats.push_back(PullStudyHelpers::meanAndStdDev(aPull));
    // stats.push_back(PullStudyHelpers::meanAndStdDev(bPull));
    // stats.push_back(PullStudyHelpers::meanAndStdDev(cPull));

    stats.push_back(PullStudyHelpers::meanAndStdDev(rPull));
    stats.push_back(PullStudyHelpers::meanAndStdDev(imZPull));

    for (auto pair = stats.begin(); pair != stats.end(); ++pair) {
        std::cout << "Pull:\t" << pair->first << "+-" << pair->second << std::endl;
    }

    PullStudyHelpers::plot_parameter_distribution("r", rPull, numExperiments);
    PullStudyHelpers::plot_parameter_distribution("imZ", imZPull, numExperiments);

    PullStudyHelpers::plotHist(chiSquaredVals, 50, "MinuitChiSq");
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
