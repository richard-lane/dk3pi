#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "../pull_study/DecaySimulator.h"
#include "../pull_study/PullStudyHelpers.h"
#include "RatioCalculator.h"
#include "fitter/MinuitPolynomialFitter.h"
#include "fitter/PhysicalFitter.h"
#include "fitter/RootFitter.h"
#include "physics.h"
#include "util.h"

/*
 * Create a dataset using accept-reject and our RatioCalculator then fit it to polynomials using ROOT's builtin TGraph
 * fitter, Minuit ChiSq and Minuit Max Likelihood
 *
 * Plot the dataset with best fit lines, and print out the fit parameters
 */
void compareRootMinuit(void)
{
    // Create an accept-reject dataset
    DecayParams_t phaseSpaceParams = {
        .x     = WORLD_AVERAGE_X,
        .y     = WORLD_AVERAGE_Y,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = 2439.0,
    };
    double maxTime = 0.002;

    size_t numCfEvents = 1e7;
    double numDcsEvents =
        PullStudyHelpers::numDCSDecays(numCfEvents, phaseSpaceParams, maxTime, 1 / phaseSpaceParams.width);

    double efficiencyTimescale = 1 / phaseSpaceParams.width;
    auto   cfRate              = [&](double x) { return Phys::cfRate(x, phaseSpaceParams, efficiencyTimescale); };
    auto   dcsRate             = [&](double x) { return Phys::dcsRate(x, phaseSpaceParams, efficiencyTimescale); };

    // Generator and PDF for random numbers
    std::random_device                     rd;
    std::shared_ptr<std::mt19937>          _gen = std::make_shared<std::mt19937>(rd());
    std::uniform_real_distribution<double> uniform;

    auto gen = [&](void) {
        double x = uniform(*_gen);
        double z = 1 - std::exp(-1 * phaseSpaceParams.width * maxTime);
        return (-1 / phaseSpaceParams.width) * std::log(1 - z * x);
    };
    auto genPDF = [&](double x) { return std::exp(-phaseSpaceParams.width * x); };

    SimulatedDecays MyDecays = SimulatedDecays(gen, genPDF, cfRate, dcsRate, std::make_pair(0., maxTime), _gen);
    MyDecays.findCfDecayTimes(numCfEvents);
    MyDecays.findDcsDecayTimes(numDcsEvents);

    // Define some time bins
    std::vector<double> timeBinLimits = util::exponentialBinLimits(maxTime, phaseSpaceParams.width, 15);

    // Divide using RatioCalculator
    std::vector<size_t> cfCounts  = util::binVector(MyDecays.RSDecayTimes, timeBinLimits);
    std::vector<size_t> dcsCounts = util::binVector(MyDecays.WSDecayTimes, timeBinLimits);
    RatioCalculator     MyRatios  = RatioCalculator(cfCounts, dcsCounts, timeBinLimits);
    MyRatios.calculateRatios();

    // Create fitters
    // Fit using ROOT builtin polynomial fit, Minuit2 fit to a + bt + ct^2
    // and Minuit2 fit to (rD2 + rD(yReZ + xImZ)Gamma t + (x2 + y2)/4 (Gamma t)2), with/without constraint on X and Y
    // All with/without integration
    FitData_t MyFitData = FitData(timeBinLimits, MyRatios.ratio, MyRatios.error);

    RootFitter BuiltInFitter = RootFitter(MyFitData);

    IntegralOptions_t      integralOptions(true, phaseSpaceParams.width, timeBinLimits, 1 / phaseSpaceParams.width);
    IntegralOptions_t      integralOptionsNoIntegral(false, 0, timeBinLimits, 0);
    MinuitPolynomialFitter MinuitPolyFit = MinuitPolynomialFitter(MyFitData, integralOptions);

    MinuitPolynomialFitter MinuitPolyFitNoIntegral = MinuitPolynomialFitter(MyFitData, integralOptionsNoIntegral);

    PhysicalFitter PhysFitter           = PhysicalFitter(MyFitData, integralOptions);
    PhysicalFitter PhysFitterNoIntegral = PhysicalFitter(MyFitData, integralOptionsNoIntegral);

    PhysicalFitter ConstrainXY           = PhysicalFitter(MyFitData, integralOptions, true);
    PhysicalFitter ConstrainXYNoIntegral = PhysicalFitter(MyFitData, integralOptionsNoIntegral, true);

    // Perform fits
    std::vector<double> initialParameterGuess{0.02, 1.0, 100.0};
    std::vector<double> initialErrorsGuess{0.01, 1.0, 100.0};
    BuiltInFitter.fit(0, maxTime, "Q");

    MinuitPolyFitNoIntegral.setPolynomialParams(initialParameterGuess, initialErrorsGuess);
    MinuitPolyFitNoIntegral.fit();

    MinuitPolyFit.setPolynomialParams(initialParameterGuess, initialErrorsGuess);
    MinuitPolyFit.fit();

    std::vector<double> initialParamGuess{phaseSpaceParams.x,
                                          phaseSpaceParams.y,
                                          phaseSpaceParams.r,
                                          phaseSpaceParams.z_im,
                                          phaseSpaceParams.z_re,
                                          phaseSpaceParams.width};
    std::vector<double> initialErrGuess{1, 1, 1, 1, 1, 1};

    // Perform a fit, fixing x, y and the width to their model values
    PhysFitter.setFitParams(initialParamGuess, initialErrGuess);
    PhysFitter.fixParameters(std::vector<std::string>{"x", "y", "width"});
    PhysFitter.fit();

    PhysFitterNoIntegral.setFitParams(initialParamGuess, initialErrGuess);
    PhysFitterNoIntegral.fixParameters(std::vector<std::string>{"x", "y", "width"});
    PhysFitterNoIntegral.fit();

    ConstrainXY.setFitParams(initialParamGuess, initialErrGuess);
    ConstrainXY.fixParameters(std::vector<std::string>{"width"});
    ConstrainXY.fit();

    ConstrainXYNoIntegral.setFitParams(initialParamGuess, initialErrGuess);
    ConstrainXYNoIntegral.fixParameters(std::vector<std::string>{"width"});
    ConstrainXYNoIntegral.fit();

    // Plot fits to file
    BuiltInFitter.plot->SetTitle("Compare Minuit and ROOT fitters;time/ns;DCS/CF ratio");
    BuiltInFitter.plot->SetLineColor(kBlack);

    MinuitPolyFit.bestFitFunction->SetLineColor(kBlue);

    MinuitPolyFitNoIntegral.bestFitFunction->SetLineColor(kTeal);
    MinuitPolyFitNoIntegral.bestFitFunction->SetLineStyle(kDotted);
    MinuitPolyFitNoIntegral.bestFitFunction->SetLineWidth(3);

    PhysFitterNoIntegral.bestFitFunction->SetLineColor(kOrange);
    PhysFitterNoIntegral.bestFitFunction->SetLineStyle(kDotted);

    PhysFitter.bestFitFunction->SetLineColor(kGreen);
    PhysFitter.bestFitFunction->SetLineStyle(kDashed);

    ConstrainXY.bestFitFunction->SetLineColor(kMagenta);
    ConstrainXY.bestFitFunction->SetLineStyle(kDashed);

    ConstrainXYNoIntegral.bestFitFunction->SetLineColor(kPink);
    ConstrainXYNoIntegral.bestFitFunction->SetLineStyle(kDashed);

    TF1*                trueFit        = new TF1("true fit", "[0] +[1]*x+[2]*x*x", 0, maxTime);
    std::vector<double> expectedParams = util::expectedParams(phaseSpaceParams);
    trueFit->SetParameter(0, expectedParams[0]);
    trueFit->SetParameter(1, expectedParams[1]);
    trueFit->SetParameter(2, expectedParams[2]);
    trueFit->SetLineColor(kGray);

    const util::LegendParams_t legendParams = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = "Compare fitters"};
    const std::vector<std::string> legendLabels{"Root best fit (red)",
                                                "Minuit polynomial fit with integral",
                                                "Minuit polynomial fit (no integral)",
                                                "Fixed x, y",
                                                "Constrained x, y",
                                                "Fixed x, y (no integral)",
                                                "Constrained x, y (no integral)",
                                                "true fit"};

    util::saveObjectsToFile<TGraph>(
        std::vector<TObject*>{BuiltInFitter.plot.get(),
                              MinuitPolyFit.bestFitFunction.get(),
                              MinuitPolyFitNoIntegral.bestFitFunction.get(),
                              PhysFitter.bestFitFunction.get(),
                              ConstrainXY.bestFitFunction.get(),
                              PhysFitterNoIntegral.bestFitFunction.get(),
                              ConstrainXYNoIntegral.bestFitFunction.get(),
                              trueFit},
        std::vector<std::string>{"AP", "CSAME", "CSAME", "CSAME", "CSAME", "CSAME", "CSAME", "CSAME"},
        legendLabels,
        "compareMinuitRootPlots.pdf",
        legendParams);

    delete trueFit;
}

int main()
{
    compareRootMinuit();
    return 0;
}
