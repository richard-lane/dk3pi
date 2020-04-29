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
#include "TH2.h"

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
 * 2d scan of im(Z) and re(Z), both with and without a constraint
 */
void test_z_scan()
{
    // Create an accept-reject dataset
    DecayParams_t phaseSpaceParams = {
        .x     = 0.0039,
        .y     = 0.0065,
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

    // Create fitters
    IntegralOptions_t integralOptions(phaseSpaceParams.width, timeBinLimits);
    FitData_t         MyFitData  = FitData(MyRatios.binCentres, MyRatios.binWidths, MyRatios.ratio, MyRatios.error);
    PhysicalFitter    PhysFitter = PhysicalFitter(MyFitData, integralOptions);
    PhysicalFitter    PhysFitterConstraint = PhysicalFitter(MyFitData, integralOptions, true);

    // Perform fit, output minimum statistic
    std::vector<double> initialParameterGuess{phaseSpaceParams.x,
                                              phaseSpaceParams.y,
                                              phaseSpaceParams.r,
                                              phaseSpaceParams.z_im,
                                              phaseSpaceParams.z_re,
                                              phaseSpaceParams.width};
    std::vector<double> initialErrorsGuess{1, 1, 1, 1, 1, 1};
    PhysFitter.setPhysicalFitParams(initialParameterGuess, initialErrorsGuess);
    PhysFitterConstraint.setPhysicalFitParams(initialParameterGuess, initialErrorsGuess);

    // Perform a 2d chi squared scan on the components of Z
    size_t numPoints      = 100;
    size_t numTotalPoints = numPoints * numPoints;
    double min            = -1;
    double max            = 1;

    PhysFitter.fixParameters(std::vector<std::string>{"width", "x", "y"});
    PhysFitter.twoDParamScan(3, 4, numPoints, numPoints, min, max, min, max);

    PhysFitterConstraint.fixParameters(std::vector<std::string>{"width"});
    PhysFitterConstraint.twoDParamScan(3, 4, numPoints, numPoints, min, max, min, max);

    std::vector<double> imVals(numTotalPoints);
    std::vector<double> reVals(numTotalPoints);
    std::vector<double> chiSquaredValsNoConstraint(numTotalPoints);
    std::vector<double> chiSquaredValsWithConstraint(numTotalPoints);

    for (size_t i = 0; i < numTotalPoints; ++i) {
        imVals[i]                       = PhysFitter.twoDParameterScan[i][0];
        reVals[i]                       = PhysFitter.twoDParameterScan[i][1];
        chiSquaredValsNoConstraint[i]   = PhysFitter.twoDParameterScan[i][2];
        chiSquaredValsWithConstraint[i] = PhysFitterConstraint.twoDParameterScan[i][2];
    }

    // Subtract off the minimum chi squared
    double minChiSqNoConstraint =
        *std::min_element(chiSquaredValsNoConstraint.begin(), chiSquaredValsNoConstraint.end());
    double minChiSqWithConstraint =
        *std::min_element(chiSquaredValsWithConstraint.begin(), chiSquaredValsWithConstraint.end());
    std::transform(chiSquaredValsNoConstraint.begin(),
                   chiSquaredValsNoConstraint.end(),
                   chiSquaredValsNoConstraint.begin(),
                   [&](double x) { return x - minChiSqNoConstraint; });
    std::transform(chiSquaredValsWithConstraint.begin(),
                   chiSquaredValsWithConstraint.end(),
                   chiSquaredValsWithConstraint.begin(),
                   [&](double x) { return x - minChiSqWithConstraint; });

    TGraph2D* NoConstraintGraph =
        new TGraph2D(numTotalPoints, imVals.data(), reVals.data(), chiSquaredValsNoConstraint.data());
    NoConstraintGraph->SetMaximum(25);
    NoConstraintGraph->SetTitle("2d Z scan;Im(Z);Re(Z);chiSq");
    util::saveObjectToFile(NoConstraintGraph, "z_.pdf", "CONT4Z");
    delete NoConstraintGraph;

    TGraph2D* ConstraintGraph =
        new TGraph2D(numTotalPoints, imVals.data(), reVals.data(), chiSquaredValsWithConstraint.data());
    ConstraintGraph->SetMaximum(25);
    ConstraintGraph->GetHistogram()->SetContour(50);
    ConstraintGraph->SetTitle("2d Z scan with constraint;Im(Z);Re(Z);chiSq");
    util::saveObjectToFile(ConstraintGraph, "z_constraint.pdf", "CONT4Z");
    delete ConstraintGraph;
}

int main()
{
    // a b c polynomial scans
    test_param_scan();
    test_2d_scan();

    // 2d scans
    test_z_scan();
    return 0;
}
