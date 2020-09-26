/* Lots of reused code here */
#include <iostream>
#include <memory>
#include <vector>

#include "DecaySimulator.h"
#include "RatioCalculator.h"
#include "fitter/MinuitPolynomialFitter.h"
#include "fitter/PhysicalFitter.h"
#include "physics.h"
#include "util.h"

#include "TCanvas.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH2.h"

#include "Minuit2/MnMinos.h"

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

    std::vector<std::pair<double, double>> parameterScan =
        MinuitChiSqFitter.chiSqParameterScan(paramIndex, numPoints, min, max);
    std::vector<double> values{};
    std::vector<double> chiSq{};
    splitVectorOfPairs(parameterScan, values, chiSq);
    TGraph* Graph = new TGraph(numPoints, values.data(), chiSq.data());
    Graph->SetTitle((graphTitle + ";value;chiSq").c_str());
    util::saveObjectToFile(Graph, (graphTitle + ".pdf").c_str(), "AP");
}

SimulatedDecays generateDecays(const FitterUtil::DecayParams_t& phaseSpaceParams,
                               const double                     maxTime,
                               const size_t                     numCfEvents,
                               const double                     efficiencyTimescale)
{
    double numDcsEvents = Phys::numDCSDecays(numCfEvents, phaseSpaceParams, maxTime, 1 / phaseSpaceParams.width);

    auto cfRate  = [&](double x) { return Phys::cfRate(x, phaseSpaceParams, efficiencyTimescale); };
    auto dcsRate = [&](double x) { return Phys::dcsRate(x, phaseSpaceParams, efficiencyTimescale); };

    // Generator and PDF for random numbers
    std::random_device                     rd;
    std::shared_ptr<std::mt19937>          _gen = std::make_shared<std::mt19937>(rd());
    std::uniform_real_distribution<double> uniform;

    auto gen = [&](void) {
        double x = uniform(*_gen);
        double z = 1 - std::exp(-1 * phaseSpaceParams.width * maxTime);
        return (-1 / phaseSpaceParams.width) * std::log(1 - z * x);
    };
    auto genPDF = [&](double x) {
        return std::exp(-phaseSpaceParams.width * x) * phaseSpaceParams.width /
               (1 - std::exp(-phaseSpaceParams.width * maxTime));
    };

    SimulatedDecays MyDecays = SimulatedDecays(gen, genPDF, cfRate, dcsRate, std::make_pair(0., maxTime), _gen);
    MyDecays.findDcsDecayTimes((size_t)numDcsEvents);
    MyDecays.findCfDecayTimes(numCfEvents);

    return MyDecays;
}

/*
 * Create a dataset using accept-reject and our RatioCalculator then fit it to polynomials using
 * Minuit ChiSq
 *
 * Perform a scan of each parameter, store them as std::vector<std::pair<double, double>> and plot
 */
void test_param_scan(void)
{
    double                    maxTime             = 0.004;
    size_t                    numCfEvents         = 700000;
    double                    efficiencyTimescale = 2500.0;
    FitterUtil::DecayParams_t phaseSpaceParams    = {
        .x     = WORLD_AVERAGE_X,
        .y     = WORLD_AVERAGE_Y,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = efficiencyTimescale,
    };

    auto MyDecays = generateDecays(phaseSpaceParams, maxTime, numCfEvents, efficiencyTimescale);

    // Define some time bins
    std::vector<double> dcsTimes{MyDecays.WSDecayTimes};
    std::sort(dcsTimes.begin(), dcsTimes.end());
    std::vector<double> timeBinLimits = util::findBinLimits(dcsTimes, 100, 0, maxTime);

    // Divide using RatioCalculator
    std::vector<size_t> cfCounts  = util::binVector(MyDecays.RSDecayTimes, timeBinLimits);
    std::vector<size_t> dcsCounts = util::binVector(MyDecays.WSDecayTimes, timeBinLimits);
    std::pair<std::vector<double>, std::vector<double>> ratiosAndErrors =
        RatioCalculator::ratioAndError(cfCounts, dcsCounts);

    // Create a fitter
    FitData_t         MyFitData = FitData(timeBinLimits, ratiosAndErrors.first, ratiosAndErrors.second);
    IntegralOptions_t integralOptions(true, phaseSpaceParams.width, timeBinLimits, efficiencyTimescale);
    PhysicalFitter    PhysicalFitter(MyFitData, integralOptions, true);

    // Perform fit, outputu minimum statistic
    std::vector<double> initialParameterGuess{phaseSpaceParams.x,
                                              phaseSpaceParams.y,
                                              phaseSpaceParams.r,
                                              phaseSpaceParams.z_im,
                                              phaseSpaceParams.z_re,
                                              phaseSpaceParams.width};
    std::vector<double> initialErrorsGuess{1, 1, 1, 1, 1, 1};
    PhysicalFitter.setFitParams(initialParameterGuess, initialErrorsGuess);
    PhysicalFitter.fixParameters(std::vector<std::string>{"width", "z_im"});
    PhysicalFitter.fit();
    std::cout << "Min chisq: " << PhysicalFitter.fitParams.fitStatistic << std::endl;
    double rSigma = PhysicalFitter.fitParams.fitParamErrors[4];
    double r0     = PhysicalFitter.fitParams.fitParams[4];

    // Minos error analysis
    size_t sigmas{3};
    PhysicalFitter._fitFcn->setErrorDef(std::sqrt((double)sigmas * (double)sigmas));
    ROOT::Minuit2::MnMinos    minos(*(PhysicalFitter._fitFcn), *(PhysicalFitter.min));
    std::pair<double, double> e0 = minos(0);
    std::cout << sigmas << "-sigma minos errors: " << std::endl;
    std::cout << "r: " << PhysicalFitter.min->UserState().Value("r") << " " << e0.first << " " << e0.second
              << std::endl;

    // Perform a chi squared scan r
    size_t                    numPoints = 400;
    std::pair<double, double> rVals(
        PhysicalFitter.fitParams.fitParams[4] - 3 * PhysicalFitter.fitParams.fitParamErrors[4],
        PhysicalFitter.fitParams.fitParams[4] + 3 * PhysicalFitter.fitParams.fitParamErrors[4]);

    std::vector<std::pair<double, double>> rScan =
        PhysicalFitter.chiSqParameterScan(4, numPoints, rVals.first, rVals.second);

    std::vector<double> values{};
    std::vector<double> chiSq{};
    splitVectorOfPairs(rScan, values, chiSq);

    // Subtract minimum chi sq
    double minChiSq = *std::min_element(chiSq.begin(), chiSq.end());
    std::transform(chiSq.begin(), chiSq.end(), chiSq.begin(), [&](double x) { return x - minChiSq; });

    TGraph* Graph = new TGraph(numPoints, values.data(), chiSq.data());
    Graph->SetTitle("3-sigma Chi2 Scan;Re(Z) value;chiSq");

    // Expected chi squared parabola
    TF1* expectedChiSq = new TF1("expected", "(x-[0])*(x-[0])/([1]*[1])", r0 - 3 * rSigma, r0 + 3 * rSigma);
    expectedChiSq->SetParameter(0, r0);
    expectedChiSq->SetParameter(1, rSigma);

    const util::LegendParams_t legend = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.9, .y2 = 0.7, .header = ""};
    util::saveObjectsToFile<TGraph>(std::vector<TObject*>{Graph, expectedChiSq},
                                    std::vector<std::string>{"AP", "CSAME"},
                                    std::vector<std::string>{"Chisq scan", "symmetric errors"},
                                    "rezScan.pdf",
                                    legend);
    delete Graph;
    delete expectedChiSq;
}

/*
 * 2d scan of im(Z) and re(Z), both with and without a constraint
 */
void test_z_scan()
{
    double                    maxTime             = 0.004;
    size_t                    numCfEvents         = 700000;
    double                    efficiencyTimescale = 2500.0;
    FitterUtil::DecayParams_t phaseSpaceParams    = {
        .x     = WORLD_AVERAGE_X,
        .y     = WORLD_AVERAGE_Y,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = efficiencyTimescale,
    };

    auto MyDecays = generateDecays(phaseSpaceParams, maxTime, numCfEvents, efficiencyTimescale);

    // Define some time bins
    std::vector<double> dcsTimes{MyDecays.WSDecayTimes};
    std::sort(dcsTimes.begin(), dcsTimes.end());
    std::vector<double> timeBinLimits = util::findBinLimits(dcsTimes, 100, 0, 1.05 * maxTime);

    // Divide using RatioCalculator
    std::vector<size_t> cfCounts  = util::binVector(MyDecays.RSDecayTimes, timeBinLimits);
    std::vector<size_t> dcsCounts = util::binVector(MyDecays.WSDecayTimes, timeBinLimits);
    std::pair<std::vector<double>, std::vector<double>> ratiosAndErrors =
        RatioCalculator::ratioAndError(cfCounts, dcsCounts);

    // Create fitters
    IntegralOptions_t integralOptions(true, phaseSpaceParams.width, timeBinLimits, 1 / phaseSpaceParams.width);
    FitData_t         MyFitData            = FitData(timeBinLimits, ratiosAndErrors.first, ratiosAndErrors.second);
    PhysicalFitter    PhysFitter           = PhysicalFitter(MyFitData, integralOptions);
    PhysicalFitter    PhysFitterConstraint = PhysicalFitter(MyFitData, integralOptions, true);

    // Perform fit, output minimum statistic
    std::vector<double> initialParameterGuess{phaseSpaceParams.x,
                                              phaseSpaceParams.y,
                                              phaseSpaceParams.r,
                                              phaseSpaceParams.z_im,
                                              phaseSpaceParams.z_re,
                                              phaseSpaceParams.width};
    std::vector<double> initialErrorsGuess{1, 1, 1, 1, 1, 1};
    PhysFitter.setFitParams(initialParameterGuess, initialErrorsGuess);
    PhysFitterConstraint.setFitParams(initialParameterGuess, initialErrorsGuess);

    // Perform a 2d chi squared scan on the components of Z
    size_t numPoints      = 100;
    size_t numTotalPoints = numPoints * numPoints;
    double min            = -1;
    double max            = 1;

    PhysFitterConstraint.fixParameters(std::vector<std::string>{"width"});
    std::vector<std::vector<double>> twoDParameterScanConstraint = PhysFitterConstraint.twoDParamScan(
        3, 4, numPoints, numPoints, std::make_pair(min, max), std::make_pair(min, max));

    std::vector<double> imVals(numTotalPoints);
    std::vector<double> reVals(numTotalPoints);
    std::vector<double> chiSquaredValsWithConstraint(numTotalPoints);

    for (size_t i = 0; i < numTotalPoints; ++i) {
        imVals[i]                       = twoDParameterScanConstraint[i][0];
        reVals[i]                       = twoDParameterScanConstraint[i][1];
        chiSquaredValsWithConstraint[i] = twoDParameterScanConstraint[i][2];
    }

    // Subtract off the minimum chi squared
    double minChiSqWithConstraint =
        *std::min_element(chiSquaredValsWithConstraint.begin(), chiSquaredValsWithConstraint.end());
    std::transform(chiSquaredValsWithConstraint.begin(),
                   chiSquaredValsWithConstraint.end(),
                   chiSquaredValsWithConstraint.begin(),
                   [&](double x) { return x - minChiSqWithConstraint; });

    // Square canvas
    TCanvas* c = new TCanvas();
    c->cd();
    c->SetWindowSize(1000, 1000);

    // Draw a circle to show the allowed Z values
    double radius = 0.56; // No idea how TEllipse works but this is the factor i need to scale my ellipse by to give it
                          // a radius of 1?
    TEllipse* boundary = new TEllipse(0, 0, radius, radius);
    boundary->SetFillColorAlpha(0, 0); // Transparent

    // Take square root to find standard deviations
    std::vector<double> constrainedFitSigmaVals{chiSquaredValsWithConstraint};
    std::transform(constrainedFitSigmaVals.begin(),
                   constrainedFitSigmaVals.end(),
                   constrainedFitSigmaVals.begin(),
                   [&](double x) { return std::sqrt(x); }); // Could use std::bind

    // Draw contour plot up to 5 sigma
    size_t    numContours                = 6;
    size_t    maxSigma                   = numContours - 1;
    double    contourLevels[numContours] = {0, 1, 2, 3, 4, 5};
    TGraph2D* ConstraintGraph =
        new TGraph2D(numTotalPoints, imVals.data(), reVals.data(), constrainedFitSigmaVals.data());
    ConstraintGraph->SetMaximum(maxSigma);
    ConstraintGraph->GetHistogram()->SetContour(numContours, contourLevels);
    ConstraintGraph->SetTitle("Z Scan;Im(Z);Re(Z)");
    ConstraintGraph->GetZaxis()->SetTitle("\\sigma");

    // Point representing "true" value of Z
    TEllipse* trueZ =
        new TEllipse(radius * phaseSpaceParams.z_im, // Again have to scale by this mysterious radius factor
                     radius * phaseSpaceParams.z_re,
                     0.01,
                     0.01);

    // Save as an image
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    ConstraintGraph->Draw("CONT4Z");
    boundary->Draw();
    trueZ->Draw();
    c->SaveAs("z_constraint.png");
    delete trueZ;
    delete boundary;
    delete c;
    delete ConstraintGraph;

    // Also plot contours in the R-delta plane
    // pretty much just a copy-paste of the above code with "magphase" stuck on the end of some variable names
    std::vector<std::pair<double, double>> magAndPhase = util::reIm2magPhase(reVals, imVals);
    std::vector<double>                    magnitudes(numTotalPoints);
    std::vector<double>                    phases(numTotalPoints);
    for (size_t i = 0; i < numTotalPoints; ++i) {
        magnitudes[i] = magAndPhase[i].first;
        phases[i]     = magAndPhase[i].second;
    }

    // Square canvas
    TCanvas* c2 = new TCanvas();
    c2->cd();
    c2->SetWindowSize(1000, 1000);

    // Draw contour plot up to 5 sigma
    size_t    numContoursMagPhase                = 6;
    size_t    maxSigmaMagPhase                   = numContoursMagPhase - 1;
    double    contourLevelsMagPhase[numContours] = {0, 1, 2, 3, 4, 5};
    TGraph2D* ConstraintGraphMagPhase =
        new TGraph2D(numTotalPoints, magnitudes.data(), phases.data(), constrainedFitSigmaVals.data());
    ConstraintGraphMagPhase->SetMaximum(maxSigmaMagPhase);
    ConstraintGraphMagPhase->GetHistogram()->SetContour(numContoursMagPhase, contourLevelsMagPhase);
    ConstraintGraphMagPhase->SetTitle("Z Scan Transformed to Polar Coords;R;\\delta;\\sigmas");

    // Point representing "true" value of Z
    std::vector<std::pair<double, double>> trueZmagPhaseVals =
        util::reIm2magPhase(std::vector<double>{phaseSpaceParams.z_re}, std::vector<double>{phaseSpaceParams.z_im});
    // Need to multiply by mysterious radius factor again
    double    trueR         = radius * (2 * trueZmagPhaseVals[0].first - 1);
    double    trueDelta     = radius * (trueZmagPhaseVals[0].second / M_PI - 1);
    TEllipse* trueZMagPhase = new TEllipse(trueR, trueDelta, 0.01, 0.01);

    // Axis limits
    ConstraintGraphMagPhase->GetXaxis()->SetLimits(0, 1);

    // Save as an image
    ConstraintGraphMagPhase->Draw("CONT4Z");
    trueZMagPhase->Draw();
    c2->SaveAs("z_Rdelta.png");
}

int main()
{
    // r scan
    test_param_scan();

    // 2d scans
    test_z_scan();
    return 0;
}