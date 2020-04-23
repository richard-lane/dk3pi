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

#include "testutil.h"

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
    std::vector<size_t> cfCounts  = util::binVector(MyDecays.RSDecayTimes, timeBinLimits);
    std::vector<size_t> dcsCounts = util::binVector(MyDecays.WSDecayTimes, timeBinLimits);
    RatioCalculator     MyRatios  = RatioCalculator(cfCounts, dcsCounts, timeBinLimits);
    MyRatios.calculateRatios();

    // Create a fitter
    FitData_t              MyFitData = FitData(MyRatios.binCentres, MyRatios.binWidths, MyRatios.ratio, MyRatios.error);
    MinuitPolynomialFitter MinuitChiSqScanner(MyFitData);

    // Perform fit, outputu minimum statistic
    std::vector<double> initialParameterGuess{0.02, 1.0, 100.0};
    std::vector<double> initialErrorsGuess{0.01, 1.0, 100.0};
    MinuitChiSqScanner.setPolynomialParams(initialParameterGuess, initialErrorsGuess);
    MinuitChiSqScanner.fit();
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
    std::vector<size_t> cfCounts  = util::binVector(MyDecays.RSDecayTimes, timeBinLimits);
    std::vector<size_t> dcsCounts = util::binVector(MyDecays.WSDecayTimes, timeBinLimits);
    RatioCalculator     MyRatios  = RatioCalculator(cfCounts, dcsCounts, timeBinLimits);
    MyRatios.calculateRatios();

    // Create a fitter
    FitData_t              MyFitData = FitData(MyRatios.binCentres, MyRatios.binWidths, MyRatios.ratio, MyRatios.error);
    MinuitPolynomialFitter MinuitChiSqScanner(MyFitData);

    // Perform fit, output minimum statistic
    std::vector<double> initialParameterGuess{0.02, 1.0, 100.0};
    std::vector<double> initialErrorsGuess{0.01, 1.0, 100.0};
    MinuitChiSqScanner.setPolynomialParams(initialParameterGuess, initialErrorsGuess);
    MinuitChiSqScanner.fit();
    std::cout << "Min chisq: " << *(MinuitChiSqScanner.statistic) << std::endl;
    std::cout << "Params: " << MinuitChiSqScanner.fitParams.fitParams[0] << " "
              << MinuitChiSqScanner.fitParams.fitParams[1] << " " << MinuitChiSqScanner.fitParams.fitParams[2]
              << std::endl;

    // Perform a 2d chi squared scan on the parameters a and b
    size_t numAPoints     = 100;
    size_t numBPoints     = 100;
    size_t numTotalPoints = numAPoints * numBPoints;
    double bSigma         = 2;
    double aSigma         = 2;
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
 * 2d scan of im(Z) and re(Z)
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
    std::vector<size_t> cfCounts  = util::binVector(MyDecays.RSDecayTimes, timeBinLimits);
    std::vector<size_t> dcsCounts = util::binVector(MyDecays.WSDecayTimes, timeBinLimits);
    RatioCalculator     MyRatios  = RatioCalculator(cfCounts, dcsCounts, timeBinLimits);
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
    PhysFitter.setPhysicalFitParams(initialParameterGuess, initialErrorsGuess);

    // Perform a 2d chi squared scan on the components of Z
    size_t numPoints      = 100;
    size_t numTotalPoints = numPoints * numPoints;
    double min            = -1;
    double max            = 1;

    PhysFitter.fixParameters(std::vector<std::string>{"width", "x", "y"});
    PhysFitter.twoDParamScan(3, 4, numPoints, numPoints, min, max, min, max);

    std::vector<double> imVals(numTotalPoints);
    std::vector<double> reVals(numTotalPoints);
    std::vector<double> chiSquaredVals(numTotalPoints);

    for (size_t i = 0; i < numTotalPoints; ++i) {
        imVals[i]         = PhysFitter.twoDParameterScan[i][0];
        reVals[i]         = PhysFitter.twoDParameterScan[i][1];
        chiSquaredVals[i] = PhysFitter.twoDParameterScan[i][2];
    }

    // Subtract off the minimum chi squared
    double minChiSq = *std::min_element(chiSquaredVals.begin(), chiSquaredVals.end());
    std::transform(
        chiSquaredVals.begin(), chiSquaredVals.end(), chiSquaredVals.begin(), [&](double x) { return x - minChiSq; });

    TGraph2D* Graph = new TGraph2D(numTotalPoints, imVals.data(), reVals.data(), chiSquaredVals.data());
    Graph->SetMaximum(25);
    Graph->SetTitle("2d Z scan;Im(Z);Re(Z);chiSq");
    util::saveObjectToFile(Graph, "z_.pdf", "CONT4Z");
    delete Graph;
}

/*
 * Perform a fit to ideal data ignore the result + perform a 2d scan of im(Z) and re(Z)
 */
void test_ideal_z_scan()
{
    std::cout << "---- Ideal dataset Z scan \n\n\n";

    // Create an ideal dataset
    double maxTime     = 0.002;
    double numTimeBins = 50;
    double error       = 0.00001;

    DecayParams_t phaseSpaceParams = {
        .x     = 0.0037,
        .y     = 0.0066,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = 2439.0,
    };
    std::vector<double> params = util::expectedParams(phaseSpaceParams);
    double              a      = params[0];
    double              b      = params[1];
    double              c      = params[2];

    double              timeBinWidth = maxTime / numTimeBins;
    std::vector<double> times        = std::vector<double>(numTimeBins, -1);
    std::vector<double> ratioErrors  = std::vector<double>(numTimeBins, -1);

    // Create idealised plot
    for (size_t i = 0; i < numTimeBins; ++i) {
        times[i] = i * timeBinWidth;
    }
    std::vector<double> ratios = idealRatios(times, error, a, b, c);
    for (size_t i = 0; i < numTimeBins; ++i) {
        ratioErrors[i] = ratio(a, b, c, times[i]) * error;
    }

    // Create a fitter + set params to an initial guess
    std::vector<double> initialParameterGuess{phaseSpaceParams.x,
                                              phaseSpaceParams.y,
                                              phaseSpaceParams.r,
                                              phaseSpaceParams.z_im,
                                              phaseSpaceParams.z_re,
                                              phaseSpaceParams.width};
    std::vector<double> initialErrorsGuess{1, 1, 1, 1, 1, 1};
    FitData_t      MyFitData  = FitData(times, std::vector<double>(ratios.size(), timeBinWidth), ratios, ratioErrors);
    PhysicalFitter PhysFitter = PhysicalFitter(MyFitData);
    PhysFitter.setPhysicalFitParams(initialParameterGuess, initialErrorsGuess);

    // PhysFitter.fit(std::vector<size_t>{3});
    // const util::LegendParams_t legend = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = ""};
    // PhysFitter.saveFitPlot("fit", "fit.pdf", &legend);

    // Perform a 2d chi squared scan on the components of Z
    size_t numPoints      = 10;
    size_t numTotalPoints = numPoints * numPoints;
    double min            = -1;
    double max            = 1;

    PhysFitter.fixParameters(std::vector<std::string>{"width", "x", "y"});
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
    util::saveObjectToFile(Graph, "z_ideal.pdf", "surf1");
    delete Graph;
}

void test_1d_z_scan()
{
    // Create an ideal dataset
    double maxTime     = 0.002;
    double numTimeBins = 50;
    double error       = 0.00001;

    DecayParams_t phaseSpaceParams = {
        .x     = 0.0037,
        .y     = 0.0066,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = 2439.0,
    };
    std::vector<double> params = util::expectedParams(phaseSpaceParams);
    double              a      = params[0];
    double              b      = params[1];
    double              c      = params[2];

    double              timeBinWidth = maxTime / numTimeBins;
    std::vector<double> times        = std::vector<double>(numTimeBins, -1);
    std::vector<double> ratioErrors  = std::vector<double>(numTimeBins, -1);

    // Create idealised plot
    for (size_t i = 0; i < numTimeBins; ++i) {
        times[i] = i * timeBinWidth;
    }
    std::vector<double> ratios = idealRatios(times, error, a, b, c);
    for (size_t i = 0; i < numTimeBins; ++i) {
        ratioErrors[i] = ratio(a, b, c, times[i]) * error;
    }

    // Create a fitter + set params to an initial guess
    std::vector<double> initialParameterGuess{phaseSpaceParams.x,
                                              phaseSpaceParams.y,
                                              phaseSpaceParams.r,
                                              phaseSpaceParams.z_im,
                                              phaseSpaceParams.z_re,
                                              phaseSpaceParams.width};
    std::vector<double> initialErrorsGuess{1, 1, 1, 1, 1, 1};
    FitData_t      MyFitData  = FitData(times, std::vector<double>(ratios.size(), timeBinWidth), ratios, ratioErrors);
    PhysicalFitter PhysFitter = PhysicalFitter(MyFitData);
    PhysFitter.setPhysicalFitParams(initialParameterGuess, initialErrorsGuess);

    // Perform a fit so we know where this dataset's parameters lie
    PhysFitter.fixParameters(std::vector<std::string>{"y", "x", "width"});
    PhysFitter.fit();
    // const util::LegendParams_t legend = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = ""};
    // PhysFitter.saveFitPlot("fit", "fit.pdf", &legend);

    // Perform 1d chi squared scans on the components of Z
    size_t numPoints = 100;
    double min       = -1;
    double max       = 1;
    // double imMin = PhysFitter.fitParams.fitParams[3] - 3 * PhysFitter.fitParams.fitParamErrors[3];
    // double imMax = PhysFitter.fitParams.fitParams[3] + 3 * PhysFitter.fitParams.fitParamErrors[3];
    // double reMin = PhysFitter.fitParams.fitParams[4] - 3 * PhysFitter.fitParams.fitParamErrors[4];
    // double reMax = PhysFitter.fitParams.fitParams[4] + 3 * PhysFitter.fitParams.fitParamErrors[4];

    // Need to fix 3 params- choose y and decay width as well as our component of Z
    PhysFitter.chiSqParameterScan(3, numPoints, min, max);

    std::vector<double> imZVals(numPoints);
    std::vector<double> reZVals(numPoints);
    std::vector<double> chiSquaredValsImScan(numPoints);
    for (size_t i = 0; i < numPoints; ++i) {
        imZVals[i]              = PhysFitter.parameterScan[i].first;
        chiSquaredValsImScan[i] = PhysFitter.parameterScan[i].second;
    }

    // Fix x and width
    PhysFitter.chiSqParameterScan(4, numPoints, min, max);
    std::vector<double> chiSquaredValsReScan(numPoints);
    for (size_t i = 0; i < numPoints; ++i) {
        reZVals[i]              = PhysFitter.parameterScan[i].first;
        chiSquaredValsReScan[i] = PhysFitter.parameterScan[i].second;
    }

    TGraph* ImGraph = new TGraph(numPoints, imZVals.data(), chiSquaredValsImScan.data());
    ImGraph->SetTitle("Im Z scan;Im(Z);chiSq");
    util::saveObjectToFile(ImGraph, "Im_Z_Scan.pdf", "");

    TGraph* ReGraph = new TGraph(numPoints, reZVals.data(), chiSquaredValsReScan.data());
    ReGraph->SetTitle("Re Z scan;Re(Z);chiSq");
    util::saveObjectToFile(ReGraph, "Re_Z_Scan.pdf", "");

    delete ImGraph;
    delete ReGraph;
}

void test_1d_r_scan()
{
    // Create an ideal dataset
    double maxTime     = 0.002;
    double numTimeBins = 50;
    double error       = 0.00001;

    DecayParams_t phaseSpaceParams = {
        .x     = 0.0037,
        .y     = 0.0066,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = 2439.0,
    };
    std::vector<double> params = util::expectedParams(phaseSpaceParams);
    double              a      = params[0];
    double              b      = params[1];
    double              c      = params[2];

    double              timeBinWidth = maxTime / numTimeBins;
    std::vector<double> times        = std::vector<double>(numTimeBins, -1);
    std::vector<double> ratioErrors  = std::vector<double>(numTimeBins, -1);

    // Create idealised plot
    for (size_t i = 0; i < numTimeBins; ++i) {
        times[i] = i * timeBinWidth;
    }
    std::vector<double> ratios = idealRatios(times, error, a, b, c);
    for (size_t i = 0; i < numTimeBins; ++i) {
        ratioErrors[i] = ratio(a, b, c, times[i]) * error;
    }

    // Create a fitter + set params to an initial guess
    std::vector<double> initialParameterGuess{phaseSpaceParams.x,
                                              phaseSpaceParams.y,
                                              phaseSpaceParams.r,
                                              phaseSpaceParams.z_im,
                                              phaseSpaceParams.z_re,
                                              phaseSpaceParams.width};
    std::vector<double> initialErrorsGuess{1, 1, 1, 1, 1, 1};
    FitData_t      MyFitData  = FitData(times, std::vector<double>(ratios.size(), timeBinWidth), ratios, ratioErrors);
    PhysicalFitter PhysFitter = PhysicalFitter(MyFitData);
    PhysFitter.setPhysicalFitParams(initialParameterGuess, initialErrorsGuess);
    PhysFitter.fixParameters(std::vector<std::string>{"z_im", "z_re", "width"});

    PhysFitter.fit();
    // const util::LegendParams_t legend = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = ""};
    // PhysFitter.saveFitPlot("fit", "fit.pdf", &legend);

    // Perform 1d chi squared scan on r
    size_t numPoints = 100;
    double min       = PhysFitter.fitParams.fitParams[0] - 3 * PhysFitter.fitParams.fitParamErrors[0];
    double max       = PhysFitter.fitParams.fitParams[0] + 3 * PhysFitter.fitParams.fitParamErrors[0];

    PhysFitter.chiSqParameterScan(0, numPoints, min, max);

    std::vector<double> rVals(numPoints);
    std::vector<double> chiSquaredVals(numPoints);
    for (size_t i = 0; i < numPoints; ++i) {
        rVals[i]          = PhysFitter.parameterScan[i].first;
        chiSquaredVals[i] = PhysFitter.parameterScan[i].second;
    }

    TGraph* Graph = new TGraph(numPoints, rVals.data(), chiSquaredVals.data());
    Graph->SetTitle("r scan;r;chiSq");
    util::saveObjectToFile(Graph, "r_Scan.pdf", "");
    delete Graph;
}

int main()
{
    // a b c polynomial scans
    test_param_scan();
    test_2d_scan();

    // 1d scans of physical fit parameters
    test_1d_r_scan();
    test_1d_z_scan();

    // 2d scans
    // test_ideal_z_scan();
    test_z_scan();
    return 0;
}
