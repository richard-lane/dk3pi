/* Lots of reused code here */
#include <iostream>
#include <vector>

#include "DecaySimulator.h"
#include "PullStudyHelpers.h"
#include "RatioCalculator.h"
#include "fitter/MinuitPolynomialFitter.h"
#include "fitter/PhysicalFitter.h"
#include "util.h"

#include "TGraph.h"
#include "TGraph2D.h"

void splitVectorOfPairs(std::vector<std::pair<double, double>>& pairs,
                        std::vector<double>&                    first,
                        std::vector<double>&                    second)
{
    for (auto it = std::make_move_iterator(pairs.begin()), end = std::make_move_iterator(pairs.end()); it != end;
         ++it) {
        first.push_back(std::move(it->first));
        second.push_back(std::move(it->second));
    }
}

void plotScan(MinuitPolynomialFitter& MinuitChiSqFitter,
              const size_t            numPoints,
              const size_t            paramIndex,
              double                  min,
              double                  max)
{
    std::string graphTitle{""};
    switch (paramIndex) {
    case 0: graphTitle = "Parameter a ChiSq scan"; break;
    case 1: graphTitle = "Parameter b ChiSq scan"; break;
    case 2: graphTitle = "Parameter c ChiSq scan"; break;
    default: throw;
    }

    MinuitChiSqFitter.chiSqParameterScan(paramIndex, numPoints, min, max);
    std::vector<double> values{};
    std::vector<double> chiSq{};
    splitVectorOfPairs(MinuitChiSqFitter.parameterScan, values, chiSq);
    TGraph* Graph = new TGraph(numPoints, values.data(), chiSq.data());
    Graph->SetTitle((graphTitle + ";value;chiSq").c_str());
    util::saveObjectToFile(Graph, (graphTitle + ".pdf").c_str(), "AP");
}

/*
 * Create a dataset using accept-reject and our RatioCalculator then fit it to polynomials using
 * Minuit ChiSq
 *
 * Perform a scan of each parameter, store them as std::vector<std::pair<double, double>> and plot
 */
void test_param_scan(void)
{
    // Create an accept-reject dataset
    DecayParams_t phaseSpaceParams = {
        .x     = 0.0037,
        .y     = 0.0066,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = 2439.0,
    };
    double maxTime      = 0.002;
    size_t numCfEvents  = 1000000;
    double numDcsEvents = PullStudyHelpers::numDCSDecays(numCfEvents, phaseSpaceParams, maxTime);

    SimulatedDecays MyDecays = SimulatedDecays(maxTime, phaseSpaceParams);
    MyDecays.findDcsDecayTimes((size_t)numDcsEvents);
    MyDecays.findCfDecayTimes(numCfEvents);

    // Define some time bins
    std::vector<double> dcsTimes{MyDecays.WSDecayTimes};
    std::sort(dcsTimes.begin(), dcsTimes.end());
    std::vector<double> timeBinLimits = util::findBinLimits(dcsTimes, 100, 0, 1.05 * maxTime);

    // Divide using RatioCalculator
    RatioCalculator MyRatios = RatioCalculator(MyDecays.RSDecayTimes, MyDecays.WSDecayTimes, timeBinLimits);
    MyRatios.calculateRatios();

    // Create a fitter
    FitData_t              MyFitData = FitData(MyRatios.binCentres, MyRatios.binWidths, MyRatios.ratio, MyRatios.error);
    MinuitPolynomialFitter MinuitChiSqScanner(MyFitData);

    // Perform fit, outputu minimum statistic
    std::vector<double> initialParameterGuess{0.02, 1.0, 100.0};
    std::vector<double> initialErrorsGuess{0.01, 1.0, 100.0};
    MinuitChiSqScanner.fit(initialParameterGuess, initialErrorsGuess, ChiSquared, std::vector<size_t>{});
    std::cout << "Min chisq: " << *(MinuitChiSqScanner.statistic) << std::endl;
    std::cout << "Params: " << MinuitChiSqScanner.fitParams.fitParams[0] << " "
              << MinuitChiSqScanner.fitParams.fitParams[1] << " " << MinuitChiSqScanner.fitParams.fitParams[2]
              << std::endl;

    // Perform a chi squared scan on each parameter
    size_t                    numPoints = 100;
    std::pair<double, double> aVals(
        MinuitChiSqScanner.fitParams.fitParams[0] - 2 * MinuitChiSqScanner.fitParams.fitParamErrors[0],
        MinuitChiSqScanner.fitParams.fitParams[0] + 2 * MinuitChiSqScanner.fitParams.fitParamErrors[0]);
    std::pair<double, double> bVals(
        MinuitChiSqScanner.fitParams.fitParams[1] - 2 * MinuitChiSqScanner.fitParams.fitParamErrors[1],
        MinuitChiSqScanner.fitParams.fitParams[1] + 2 * MinuitChiSqScanner.fitParams.fitParamErrors[1]);
    std::pair<double, double> cVals(
        MinuitChiSqScanner.fitParams.fitParams[2] - 2 * MinuitChiSqScanner.fitParams.fitParamErrors[2],
        MinuitChiSqScanner.fitParams.fitParams[2] + 2 * MinuitChiSqScanner.fitParams.fitParamErrors[2]);

    plotScan(MinuitChiSqScanner, numPoints, 0, aVals.first, aVals.second);
    plotScan(MinuitChiSqScanner, numPoints, 1, bVals.first, bVals.second);
    plotScan(MinuitChiSqScanner, numPoints, 2, cVals.first, cVals.second);
}

void test_2d_scan()
{
    // Create an accept-reject dataset
    DecayParams_t phaseSpaceParams = {
        .x     = 0.0037,
        .y     = 0.0066,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = 2439.0,
    };
    double maxTime = 0.002;

    size_t numCfEvents  = 1000000;
    double numDcsEvents = PullStudyHelpers::numDCSDecays(numCfEvents, phaseSpaceParams, maxTime);

    SimulatedDecays MyDecays = SimulatedDecays(maxTime, phaseSpaceParams);
    MyDecays.findDcsDecayTimes((size_t)numDcsEvents);
    MyDecays.findCfDecayTimes(numCfEvents);

    // Define some time bins
    std::vector<double> dcsTimes{MyDecays.WSDecayTimes};
    std::sort(dcsTimes.begin(), dcsTimes.end());
    std::vector<double> timeBinLimits = util::findBinLimits(dcsTimes, 100, 0, 1.05 * maxTime);

    // Divide using RatioCalculator
    RatioCalculator MyRatios = RatioCalculator(MyDecays.RSDecayTimes, MyDecays.WSDecayTimes, timeBinLimits);
    MyRatios.calculateRatios();

    // Create a fitter
    FitData_t              MyFitData = FitData(MyRatios.binCentres, MyRatios.binWidths, MyRatios.ratio, MyRatios.error);
    MinuitPolynomialFitter MinuitChiSqScanner(MyFitData);

    // Perform fit, output minimum statistic
    std::vector<double> initialParameterGuess{0.02, 1.0, 100.0};
    std::vector<double> initialErrorsGuess{0.01, 1.0, 100.0};
    MinuitChiSqScanner.fit(initialParameterGuess, initialErrorsGuess, ChiSquared, std::vector<size_t>{});
    std::cout << "Min chisq: " << *(MinuitChiSqScanner.statistic) << std::endl;
    std::cout << "Params: " << MinuitChiSqScanner.fitParams.fitParams[0] << " "
              << MinuitChiSqScanner.fitParams.fitParams[1] << " " << MinuitChiSqScanner.fitParams.fitParams[2]
              << std::endl;

    // Perform a 2d chi squared scan on the parameters a and b
    size_t numAPoints     = 100;
    size_t numBPoints     = 100;
    size_t numTotalPoints = numAPoints * numBPoints;
    double bSigma         = 3;
    double aSigma         = 10;
    double minA = MinuitChiSqScanner.fitParams.fitParams[0] - aSigma * MinuitChiSqScanner.fitParams.fitParamErrors[0];
    double minB = MinuitChiSqScanner.fitParams.fitParams[1] - bSigma * MinuitChiSqScanner.fitParams.fitParamErrors[1];
    double maxA = MinuitChiSqScanner.fitParams.fitParams[0] + aSigma * MinuitChiSqScanner.fitParams.fitParamErrors[0];
    double maxB = MinuitChiSqScanner.fitParams.fitParams[1] + bSigma * MinuitChiSqScanner.fitParams.fitParamErrors[1];

    MinuitChiSqScanner.twoDParamScan(0, 1, numAPoints, numBPoints, minA, maxA, minB, maxB);

    std::vector<double> aVals(numTotalPoints);
    std::vector<double> bVals(numTotalPoints);
    std::vector<double> chiSquaredVals(numTotalPoints);

    for (size_t i = 0; i < numTotalPoints; ++i) {
        aVals[i]          = MinuitChiSqScanner.twoDParameterScan[i][0];
        bVals[i]          = MinuitChiSqScanner.twoDParameterScan[i][1];
        chiSquaredVals[i] = MinuitChiSqScanner.twoDParameterScan[i][2];
    }

    TGraph2D* Graph = new TGraph2D(numTotalPoints, aVals.data(), bVals.data(), chiSquaredVals.data());
    Graph->SetTitle("2d a b scan;a;b;chiSq");
    util::saveObjectToFile(Graph, "a_.pdf", "surf1");
    delete Graph;
    std::cout << std::endl;
}

/*
 * Perform a fit, ignore the result + perform a 2d scan of im(Z) and re(Z)
 */
void test_z_scan()
{
    // Create an accept-reject dataset
    DecayParams_t phaseSpaceParams = {
        .x     = 0.0037,
        .y     = 0.0066,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = 2439.0,
    };
    double maxTime = 0.002;

    size_t numCfEvents  = 1000000;
    double numDcsEvents = PullStudyHelpers::numDCSDecays(numCfEvents, phaseSpaceParams, maxTime);

    SimulatedDecays MyDecays = SimulatedDecays(maxTime, phaseSpaceParams);
    MyDecays.findDcsDecayTimes((size_t)numDcsEvents);
    MyDecays.findCfDecayTimes(numCfEvents);

    // Define some time bins
    std::vector<double> dcsTimes{MyDecays.WSDecayTimes};
    std::sort(dcsTimes.begin(), dcsTimes.end());
    std::vector<double> timeBinLimits = util::findBinLimits(dcsTimes, 100, 0, 1.05 * maxTime);

    // Divide using RatioCalculator
    RatioCalculator MyRatios = RatioCalculator(MyDecays.RSDecayTimes, MyDecays.WSDecayTimes, timeBinLimits);
    MyRatios.calculateRatios();

    // Create a fitter
    FitData_t      MyFitData  = FitData(MyRatios.binCentres, MyRatios.binWidths, MyRatios.ratio, MyRatios.error);
    PhysicalFitter PhysFitter = PhysicalFitter(MyFitData);

    // Perform fit, output minimum statistic
    std::vector<double> initialParameterGuess{phaseSpaceParams.x,
                                              phaseSpaceParams.y,
                                              phaseSpaceParams.r,
                                              phaseSpaceParams.z_im,
                                              phaseSpaceParams.z_re,
                                              phaseSpaceParams.width};
    std::vector<double> initialErrorsGuess{1, 1, 1, 1, 1, 1};
    PhysFitter.fit(initialParameterGuess, initialErrorsGuess, ChiSquared, std::vector<size_t>{0});
    std::cout << "Min chisq: " << *(PhysFitter.statistic) << std::endl;

    // Perform a 2d chi squared scan on the components of Z
    size_t numPoints      = 100;
    size_t numTotalPoints = numPoints * numPoints;
    double min            = -1;
    double max            = 1;

    PhysFitter.twoDParamScan(3, 4, numPoints, numPoints, min, max, min, max);

    std::vector<double> imVals(numTotalPoints);
    std::vector<double> reVals(numTotalPoints);
    std::vector<double> chiSquaredVals(numTotalPoints);

    for (size_t i = 0; i < numTotalPoints; ++i) {
        imVals[i]         = PhysFitter.twoDParameterScan[i][0];
        reVals[i]         = PhysFitter.twoDParameterScan[i][1];
        chiSquaredVals[i] = PhysFitter.twoDParameterScan[i][2];
    }

    TGraph2D* Graph = new TGraph2D(numTotalPoints, imVals.data(), reVals.data(), chiSquaredVals.data());
    Graph->SetTitle("2d Z scan;Im(Z);Re(Z);chiSq");
    util::saveObjectToFile(Graph, "z_.pdf", "surf1");
    delete Graph;
}

int main()
{
    test_param_scan();
    test_2d_scan();
    test_z_scan();
    return 0;
}
