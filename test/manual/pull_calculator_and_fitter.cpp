/*
 * Input ideal decay counts our calculator and fitter and hope that it outputs the right fit parameters
 */

#include <iostream>
#include <random>
#include <vector>

#include "../pull_study/DecaySimulator.h"
#include "../pull_study/PullStudyHelpers.h"
#include "RatioCalculator.h"
#include "fitter/MinuitPolynomialFitter.h"
#include "physics.h"

#include "TH1I.h"

#include <boost/progress.hpp>

#include "testutil.h"

/*
 * Find the number of DCS points expected in a bin
 *
 * Must provide the integral of the DCS rate between 0 maxtime
 */
size_t dcsPointsInBin(const double         binLow,
                      const double         binHigh,
                      const size_t         totalPoints,
                      const double         dcsIntegral,
                      const DecayParams_t &decayParams)
{
    double binIntegral = Phys::dcsIntegral(binLow, binHigh, decayParams, 1e-12, 12);
    return (size_t)(totalPoints * binIntegral / dcsIntegral);
}

/*
 * Find the number of CF points expected in a bin
 *
 * Must provide the integral of the CF rate between 0 and maxtime
 */
size_t cfPointsInBin(const double         binLow,
                     const double         binHigh,
                     const size_t         totalPoints,
                     const double         cfIntegral,
                     const DecayParams_t &decayParams)
{
    double binIntegral = Phys::cfIntegral(binLow, binHigh, decayParams, 1e-12, 12);
    return (size_t)(totalPoints * binIntegral / cfIntegral);
}

/*
 * Find the number of DCS points in each bin, drawing the counts from poisson distributions
 */
std::vector<size_t> findDcsCounts(const std::vector<double> &binLimits,
                                  const size_t               totalNumEvents,
                                  const double               dcsIntegral,
                                  const DecayParams_t &      decayParams,
                                  std::mt19937 &             gen)
{
    std::vector<size_t> counts(binLimits.size() - 1);

    for (size_t i = 0; i < counts.size(); ++i) {
        double meanCounts = dcsPointsInBin(binLimits[i], binLimits[i + 1], totalNumEvents, dcsIntegral, decayParams);
        std::poisson_distribution<size_t> dist(meanCounts);
        counts[i] = dist(gen);
    }
    return counts;
}

/*
 * Find the number of CF points in each bin, drawing the counts from poisson distributions
 */
std::vector<size_t> findCfCounts(const std::vector<double> &binLimits,
                                 const size_t               totalNumEvents,
                                 const double               cfIntegral,
                                 const DecayParams_t &      decayParams,
                                 std::mt19937 &             gen)
{
    std::vector<size_t> counts(binLimits.size() - 1);

    for (size_t i = 0; i < counts.size(); ++i) {
        double meanCounts = cfPointsInBin(binLimits[i], binLimits[i + 1], totalNumEvents, cfIntegral, decayParams);
        std::poisson_distribution<size_t> dist(meanCounts);
        counts[i] = dist(gen);
    }
    return counts;
}

/*
 * Create idealised binned DCS and CF counts, divide them with a RatioCalculator + fit the ratios with a
 * MinuitPolynomialFitter
 *
 * Takes loads of parameters because i cba to make it better
 *
 * returns the fitter cus why not
 */
MinuitPolynomialFitter performFit(const DecayParams_t &      decayParams,
                                  const double               maxTime,
                                  const size_t               numDcsEvents,
                                  const size_t               numCfEvents,
                                  const std::vector<double> &binLimits)
{
    // Create a random number generator that will be used to find the number of counts in each bin
    std::random_device rd;
    std::mt19937       gen(rd());

    // Find the counts in each bin
    const double        dcsIntegral = Phys::dcsIntegral(0, maxTime, decayParams);
    const double        cfIntegral  = Phys::cfIntegral(0, maxTime, decayParams);
    std::vector<size_t> dcsCounts   = findDcsCounts(binLimits, numDcsEvents, dcsIntegral, decayParams, gen);
    std::vector<size_t> cfCounts    = findCfCounts(binLimits, numCfEvents, cfIntegral, decayParams, gen);

    // Divide
    RatioCalculator MyRatioCalculator(cfCounts, dcsCounts, binLimits);
    MyRatioCalculator.calculateRatios();

    // Fit our idealised plot
    FitData_t MyFitData(
        MyRatioCalculator.binCentres, MyRatioCalculator.binWidths, MyRatioCalculator.ratio, MyRatioCalculator.error);
    MinuitPolynomialFitter MyFitter(MyFitData);

    MyFitter.setPolynomialParams(util::expectedParams(decayParams), std::vector<double>{1, 1, 1});
    MyFitter.fit();

    return MyFitter;
}

int main()
{
    DecayParams_t DecayParams = {
        .x     = 0.004,
        .y     = 0.007,
        .r     = 0.05,
        .z_im  = -0.3,
        .z_re  = 0.8,
        .width = 2500.0,
    };

    double maxTime     = 0.005;
    size_t numTimeBins = 50;

    std::vector<double> binLimits = util::exponentialBinLimits(maxTime, DecayParams.width, numTimeBins);

    double numCfEvents  = 10e6;
    double numDcsEvents = PullStudyHelpers::numDCSDecays(numCfEvents, DecayParams, maxTime);

    size_t              numExperiments = 1000;
    std::vector<double> aPull          = std::vector<double>(numExperiments, std::nan("-1"));
    std::vector<double> bPull          = std::vector<double>(numExperiments, std::nan("-1"));
    std::vector<double> cPull          = std::vector<double>(numExperiments, std::nan("-1"));

    std::vector<double> params = util::expectedParams(DecayParams);
    double              a      = params[0];
    double              b      = params[1];
    double              c      = params[2];

    boost::progress_display showProgress(numExperiments);

    for (size_t i = 0; i < numExperiments; ++i) {
        MinuitPolynomialFitter MyFitter  = performFit(DecayParams, maxTime, numDcsEvents, numCfEvents, binLimits);
        std::vector<double>    outParams = MyFitter.fitParams.fitParams;
        std::vector<double>    outErrors = MyFitter.fitParams.fitParamErrors;
        aPull[i]                         = (outParams[0] - a) / outErrors[0];
        bPull[i]                         = (outParams[1] - b) / outErrors[1];
        cPull[i]                         = (outParams[2] - c) / outErrors[2];
        ++showProgress;

        // const util::LegendParams_t legend   = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = ""};
        // MyFitter.saveFitPlot("Ratio fit from generated counts", "test_calc_fit.pdf", &legend);
    }

    PullStudyHelpers::plot_parameter_distribution("a", aPull, numExperiments);
    PullStudyHelpers::plot_parameter_distribution("b", bPull, numExperiments);
    PullStudyHelpers::plot_parameter_distribution("c", cPull, numExperiments);

    return 0;
}
