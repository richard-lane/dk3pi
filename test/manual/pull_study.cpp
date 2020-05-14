#include <iostream>

#include "DecaySimulator.h"
#include "MinuitFcns.h"
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
        .x     = WORLD_AVERAGE_X,
        .y     = WORLD_AVERAGE_Y,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = 2439.0,
    };
    double maxTime = 10 / phaseSpaceParams.width;
    // std::vector<double> expectedFitParams = util::expectedParams(phaseSpaceParams);
    size_t numBins = 100;

    // Create RNGs for numbers of decays
    double       meanNumDcsDecays = PullStudyHelpers::numDCSDecays(meanNumCfEvents, phaseSpaceParams, maxTime);
    std::mt19937 gen;
    std::poisson_distribution<size_t> cfDist(meanNumCfEvents);
    std::poisson_distribution<size_t> dcsDist(meanNumDcsDecays);

    // Find exponentially-spaced time bin limits to use
    // std::vector<double> binLimits = util::exponentialBinLimits(maxTime, phaseSpaceParams.width, numBins);
    std::vector<double> binLimits(numBins + 1);
    for (size_t i = 0; i <= numBins; ++i) {
        binLimits[i] = i * maxTime / numBins;
    }

    // Create a decay simulator
    SimulatedDecays MyDecays(maxTime, phaseSpaceParams);

    // Initialise vectors of fit parameter pulls and chi squared
    // std::vector<double> aPull(numExperiments, -1);
    // std::vector<double> bPull(numExperiments, -1);
    // std::vector<double> cPull(numExperiments, -1);
    std::vector<double> rPull(numExperiments, -1);
    std::vector<double> xPull(numExperiments, -1);
    std::vector<double> yPull(numExperiments, -1);
    std::vector<double> chiSquaredVals(numExperiments, -1);

    boost::progress_display showProgress(numExperiments);

    for (size_t i = 0; i < numExperiments; ++i) {

        // Find how many decays to simulate
        size_t numCfEvents  = cfDist(gen);
        size_t numDcsEvents = dcsDist(gen);

        // Simulate them
        MyDecays.findCfDecayTimes(numCfEvents);
        MyDecays.findDcsDecayTimes(numDcsEvents);
        MyDecays.plotRates(binLimits);

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

        PhysicalFitter MyFitter(MyFitData, *integralOptions, false);
        MyFitter.setPhysicalFitParams(std::vector<double>{phaseSpaceParams.x,
                                                          phaseSpaceParams.y,
                                                          phaseSpaceParams.r,
                                                          phaseSpaceParams.z_im,
                                                          phaseSpaceParams.z_re,
                                                          phaseSpaceParams.width},
                                      std::vector<double>(6, 1));
        MyFitter.fixParameters(std::vector<std::string>{"width", "z_re", "z_im"});
        MyFitter.fit();
        // const util::LegendParams_t legend = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = ""};
        // MyFitter.saveFitPlot("fit", "fit.pdf", &legend);

        // Store parameter and chi squared
        xPull[i] = (MyFitter.fitParams.fitParams[0] - phaseSpaceParams.x) / MyFitter.fitParams.fitParamErrors[0];
        yPull[i] = (MyFitter.fitParams.fitParams[1] - phaseSpaceParams.y) / MyFitter.fitParams.fitParamErrors[1];
        rPull[i] = (MyFitter.fitParams.fitParams[2] - phaseSpaceParams.r) / MyFitter.fitParams.fitParamErrors[2];
        chiSquaredVals[i] = *MyFitter.statistic;

        ++showProgress;
    }

    // Output mean and std dev of pulls
    std::vector<std::pair<double, double>> stats{};
    // stats.push_back(PullStudyHelpers::meanAndStdDev(aPull));
    // stats.push_back(PullStudyHelpers::meanAndStdDev(bPull));
    // stats.push_back(PullStudyHelpers::meanAndStdDev(cPull));

    stats.push_back(PullStudyHelpers::meanAndStdDev(rPull));
    stats.push_back(PullStudyHelpers::meanAndStdDev(xPull));
    stats.push_back(PullStudyHelpers::meanAndStdDev(yPull));

    for (auto pair = stats.begin(); pair != stats.end(); ++pair) {
        std::cout << "Pull:\t" << pair->first << "+-" << pair->second << std::endl;
    }

    PullStudyHelpers::plot_parameter_distribution("r", rPull, numExperiments);
    PullStudyHelpers::plot_parameter_distribution("x", xPull, numExperiments);
    PullStudyHelpers::plot_parameter_distribution("y", yPull, numExperiments);

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

    pull_study(1e7, 100, integrate);
}
