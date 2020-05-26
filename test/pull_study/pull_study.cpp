#include <iostream>

#include "DecaySimulator.h"
#include "MinuitPolynomialFitter.h"
#include "PhysicalFitter.h"
#include "RatioCalculator.h"

#include "FitterUtils.h"
#include "PullStudyHelpers.h"
#include "util.h"

#include <boost/progress.hpp>

void pull_study(const size_t meanNumCfEvents, const size_t numExperiments)
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
    double              maxTime             = 10 / phaseSpaceParams.width;
    std::vector<double> expectedFitParams   = util::expectedParams(phaseSpaceParams);
    size_t              numBins             = 25;
    double              efficiencyTimescale = 1 / phaseSpaceParams.width;

    // Create RNGs for numbers of decays
    double meanNumDcsDecays =
        PullStudyHelpers::numDCSDecays(meanNumCfEvents, phaseSpaceParams, maxTime, efficiencyTimescale);
    std::mt19937                      gen;
    std::poisson_distribution<size_t> cfDist(meanNumCfEvents);
    std::poisson_distribution<size_t> dcsDist(meanNumDcsDecays);

    // Find exponentially-spaced time bin limits to use
    std::vector<double> binLimits = util::exponentialBinLimits(maxTime, phaseSpaceParams.width, numBins);

    // Make some bins at the start wider (because of the efficiency)
    binLimits.erase(binLimits.begin() + 1, binLimits.begin() + 3);
    binLimits.erase(binLimits.begin() + 4);
    // std::vector<double> binLimits(numBins + 1);
    // for (size_t i = 0; i <= numBins; ++i) {
    //     binLimits[i] = i * maxTime / numBins;
    // }

    // Create a decay simulator
    SimulatedDecays MyDecays(maxTime, phaseSpaceParams, efficiencyTimescale);

    // Initialise vectors of fit parameter pulls and chi squared
    std::vector<double> aPull(numExperiments, -1);
    std::vector<double> bPull(numExperiments, -1);
    std::vector<double> cPull(numExperiments, -1);
    std::vector<double> polyChiSqVals(numExperiments, -1);

    std::vector<double> rPull(numExperiments, -1);
    std::vector<double> reZPull(numExperiments, -1);
    std::vector<double> chiSquaredVals(numExperiments, -1);

    boost::progress_display showProgress(numExperiments);

    for (size_t i = 0; i < numExperiments; ++i) {

        // Find how many decays to simulate
        size_t numCfEvents  = cfDist(gen);
        size_t numDcsEvents = dcsDist(gen);

        // Simulate them
        MyDecays.findCfDecayTimes(numCfEvents);
        MyDecays.findDcsDecayTimes(numDcsEvents);
        // MyDecays.plotRates(binLimits);

        // Time binning
        std::vector<size_t> cfCounts  = util::binVector(MyDecays.RSDecayTimes, binLimits);
        std::vector<size_t> dcsCounts = util::binVector(MyDecays.WSDecayTimes, binLimits);

        // Find the ratio
        RatioCalculator MyRatios(cfCounts, dcsCounts, binLimits);
        MyRatios.calculateRatios();

        IntegralOptions_t integralOptions(true, phaseSpaceParams.width, binLimits, efficiencyTimescale);
        FitData_t         MyFitData(binLimits, MyRatios.ratio, MyRatios.error);

        // Fit data with minuit polynomial fitter and with constrained X, Y
        MinuitPolynomialFitter PolyFitter(MyFitData, integralOptions);
        PolyFitter.setPolynomialParams(expectedFitParams, std::vector<double>(3, 1));
        PolyFitter.fit();

        PhysicalFitter MyFitter(MyFitData, integralOptions, true);
        MyFitter.setFitParams(std::vector<double>{phaseSpaceParams.x,
                                                  phaseSpaceParams.y,
                                                  phaseSpaceParams.r,
                                                  phaseSpaceParams.z_im,
                                                  phaseSpaceParams.z_re,
                                                  phaseSpaceParams.width},
                              std::vector<double>(6, 1));
        MyFitter.fixParameters(std::vector<std::string>{"width", "z_im"});
        MyFitter.fit();
        // const util::LegendParams_t legend = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = ""};
        // MyFitter.saveFitPlot("fit", "fit.pdf", &legend);

        // Store parameter and chi squared
        aPull[i] = (PolyFitter.fitParams.fitParams[0] - expectedFitParams[0]) / PolyFitter.fitParams.fitParamErrors[0];
        bPull[i] = (PolyFitter.fitParams.fitParams[1] - expectedFitParams[1]) / PolyFitter.fitParams.fitParamErrors[1];
        cPull[i] = (PolyFitter.fitParams.fitParams[2] - expectedFitParams[2]) / PolyFitter.fitParams.fitParamErrors[2];
        polyChiSqVals[i] = PolyFitter.fitParams.fitStatistic;

        reZPull[i] = (MyFitter.fitParams.fitParams[4] - phaseSpaceParams.z_re) / MyFitter.fitParams.fitParamErrors[4];
        rPull[i]   = (MyFitter.fitParams.fitParams[2] - phaseSpaceParams.r) / MyFitter.fitParams.fitParamErrors[2];
        chiSquaredVals[i] = MyFitter.fitParams.fitStatistic;

        ++showProgress;
    }

    // Output mean and std dev of pulls
    std::vector<std::pair<double, double>> polyStats{};
    std::vector<std::pair<double, double>> stats{};

    polyStats.push_back(PullStudyHelpers::meanAndStdDev(aPull));
    polyStats.push_back(PullStudyHelpers::meanAndStdDev(bPull));
    polyStats.push_back(PullStudyHelpers::meanAndStdDev(cPull));

    stats.push_back(PullStudyHelpers::meanAndStdDev(rPull));
    stats.push_back(PullStudyHelpers::meanAndStdDev(reZPull));

    std::cout << "Constrained XY Fit:" << std::endl;
    for (auto pair = stats.begin(); pair != stats.end(); ++pair) {
        std::cout << "Pull:\t" << pair->first << "+-" << pair->second << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Polynomial fit:" << std::endl;
    for (auto pair = polyStats.begin(); pair != polyStats.end(); ++pair) {
        std::cout << "Pull:\t" << pair->first << "+-" << pair->second << std::endl;
    }
    std::cout << std::endl;

    PullStudyHelpers::plot_parameter_distribution("r", rPull, numExperiments);
    PullStudyHelpers::plot_parameter_distribution("Re(Z)", reZPull, numExperiments);

    PullStudyHelpers::plot_parameter_distribution("a", aPull, numExperiments);
    PullStudyHelpers::plot_parameter_distribution("b", bPull, numExperiments);
    PullStudyHelpers::plot_parameter_distribution("c", cPull, numExperiments);

    PullStudyHelpers::plotHist(chiSquaredVals, 50, "MinuitChiSq");
    PullStudyHelpers::plotHist(polyChiSqVals, 50, "polyMinuitChiSq");
}

int main()
{
    pull_study(1e6, 100);
}
