/*
 * Input ideal decay counts our calculator and fitter and hope that it outputs the right fit parameters
 */

#include <iostream>
#include <random>
#include <vector>

#include "DecaySimulator.h"
#include "PullStudyHelpers.h"
#include "RatioCalculator.h"
#include "fitter/MinuitPolynomialFitter.h"

#include "TH1I.h"

#include "testutil.h"

/*
 * Find the number of DCS points expected in a bin
 */
size_t dcsPointsInBin(const double         binLow,
                      const double         binHigh,
                      const size_t         totalPoints,
                      const double         maxTime,
                      const DecayParams_t &decayParams)
{
    double totalIntegral = util::dcsIntegral(0, maxTime, decayParams);
    double binIntegral   = util::dcsIntegral(binLow, binHigh, decayParams);
    return (size_t)(totalPoints * binIntegral / totalIntegral);
}

/*
 * Find the number of CF points expected in a bin
 */
size_t cfPointsInBin(const double         binLow,
                     const double         binHigh,
                     const size_t         totalPoints,
                     const double         maxTime,
                     const DecayParams_t &decayParams)
{
    double totalIntegral = util::cfIntegral(0, maxTime, decayParams);
    double binIntegral   = util::cfIntegral(binLow, binHigh, decayParams);
    return (size_t)(totalPoints * binIntegral / totalIntegral);
}

/*
 * Find the number of DCS points in each bin
 */
std::vector<size_t> findDcsCounts(const std::vector<double> &binLimits,
                                  const size_t               totalNumEvents,
                                  const double               maxTime,
                                  const DecayParams_t &      decayParams)
{
    std::vector<size_t> counts(binLimits.size() - 1);

    for (size_t i = 0; i < counts.size(); ++i) {
        counts[i] = dcsPointsInBin(binLimits[i], binLimits[i + 1], totalNumEvents, maxTime, decayParams);
    }
    return counts;
}

/*
 * Find the number of CF points in each bin
 */
std::vector<size_t> findCfCounts(const std::vector<double> &binLimits,
                                 const size_t               totalNumEvents,
                                 const double               maxTime,
                                 const DecayParams_t &      decayParams)
{
    std::vector<size_t> counts(binLimits.size() - 1);

    for (size_t i = 0; i < counts.size(); ++i) {
        counts[i] = cfPointsInBin(binLimits[i], binLimits[i + 1], totalNumEvents, maxTime, decayParams);
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
    // Find the counts in each bin
    std::vector<size_t> dcsCounts = findDcsCounts(binLimits, numDcsEvents, maxTime, decayParams);
    std::vector<size_t> cfCounts  = findCfCounts(binLimits, numCfEvents, maxTime, decayParams);

    // TH1I *DcsHist = new TH1I("dcs", "dcs", binLimits.size() - 1, 0, maxTime);
    // TH1I *CfHist  = new TH1I("cf", "cf", binLimits.size() - 1, 0, maxTime);
    // for (size_t i = 0; i < binLimits.size() - 1; i++) {
    //    DcsHist->SetBinContent(i + 1, dcsCounts[i]);
    //    CfHist->SetBinContent(i + 1, cfCounts[i]);
    //}
    // util::saveObjectToFile(DcsHist, "dcs.pdf", "");
    // util::saveObjectToFile(CfHist, "cf.pdf", "");
    // delete DcsHist;
    // delete CfHist;

    // Divide
    RatioCalculator MyRatioCalculator(cfCounts, dcsCounts, binLimits);
    MyRatioCalculator.calculateRatios();

    // Fit our idealised plot
    FitData_t MyFitData(
        MyRatioCalculator.binCentres, MyRatioCalculator.binWidths, MyRatioCalculator.ratio, MyRatioCalculator.error);
    MinuitPolynomialFitter MyFitter(MyFitData, binLimits, decayParams.width);

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

    std::vector<double> binLimits = PullStudyHelpers::exponentialBinLimits(maxTime, DecayParams.width, numTimeBins);

    double numCfEvents  = 10e6;
    double numDcsEvents = PullStudyHelpers::numDCSDecays(numCfEvents, DecayParams, maxTime);

    // Find idealised counts in each bin and perform a fit. It should fit Perfectly
    MinuitPolynomialFitter     MyFitter = performFit(DecayParams, maxTime, numDcsEvents, numCfEvents, binLimits);
    const util::LegendParams_t legend   = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = ""};
    MyFitter.saveFitPlot("Ratio fit from generated counts", "test_calc_fit.pdf", &legend);

    std::cout << "Fit statistic:\t" << *MyFitter.statistic << std::endl;

    return 0;
}
