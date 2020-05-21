#include <iostream>
#include <utility>
#include <vector>

#include "../pull_study/DecaySimulator.h"
#include "../pull_study/PullStudyHelpers.h"
#include "RatioCalculator.h"
#include "fitter/MinuitPolynomialFitter.h"
#include "fitter/PhysicalFitter.h"
#include "fitter/RootFitter.h"
#include "util.h"

/*
 * Create a dataset using accept-reject and our RatioCalculator then fit it to polynomials using ROOT's builtin TGraph
 * fitter, Minuit ChiSq and Minuit Max Likelihood
 *
 * Plot the dataset along with all three trendlines, and print out the fit parameters
 */
void compareRootMinuit(void)
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

    SimulatedDecays MyDecays = SimulatedDecays(maxTime, phaseSpaceParams, 0);
    MyDecays.findDcsDecayTimes((size_t)numDcsEvents);
    MyDecays.findCfDecayTimes(numCfEvents);

    // Define some time bins
    std::vector<double> dcsTimes{MyDecays.WSDecayTimes};
    std::sort(dcsTimes.begin(), dcsTimes.end());
    std::vector<double> timeBinLimits = util::exponentialBinLimits(maxTime, phaseSpaceParams.width, 10);

    // Divide using RatioCalculator
    std::vector<size_t> cfCounts  = util::binVector(MyDecays.RSDecayTimes, timeBinLimits);
    std::vector<size_t> dcsCounts = util::binVector(MyDecays.WSDecayTimes, timeBinLimits);
    RatioCalculator     MyRatios  = RatioCalculator(cfCounts, dcsCounts, timeBinLimits);
    MyRatios.calculateRatios();

    // Create 4 fitters
    FitData_t MyFitData = FitData(MyRatios.binCentres, MyRatios.binWidths, MyRatios.ratio, MyRatios.error);

    // Root builtin fitter
    RootFitter CernFitter = RootFitter(MyFitData);

    // Polynomial fit with integration
    IntegralOptions_t      integralOptions(phaseSpaceParams.width, timeBinLimits, 1e-10, 10);
    MinuitPolynomialFitter MinuitPolyFit = MinuitPolynomialFitter(MyFitData, &integralOptions);

    // Polynomial fit without integration
    MinuitPolynomialFitter MinuitPolyFitNoIntegral = MinuitPolynomialFitter(MyFitData);

    // Fit to phase space params
    PhysicalFitter PhysFitter = PhysicalFitter(MyFitData, integralOptions);

    // Perform fits
    std::vector<double> initialParameterGuess{0.02, 1.0, 100.0};
    std::vector<double> initialErrorsGuess{0.01, 1.0, 100.0};
    CernFitter.fit(0, maxTime, "Q");

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
    PhysFitter.setPhysicalFitParams(initialParamGuess, initialErrGuess);
    PhysFitter.fixParameters(std::vector<std::string>{"x", "y", "width"});
    PhysFitter.fit();

    // Print fit parameters to console
    for (size_t i = 0; i < 3; ++i) {
        std::cout << "ROOT fit params: " << CernFitter.fitParams.fitParams[i] << "+-"
                  << CernFitter.fitParams.fitParamErrors[i] << std::endl;
    }
    std::cout << "\tChiSq = " << *(CernFitter.statistic) << std::endl;

    for (size_t i = 0; i < 3; ++i) {
        std::cout << "Polyfit params: " << MinuitPolyFit.fitParams.fitParams[i] << "+-"
                  << MinuitPolyFit.fitParams.fitParamErrors[i] << std::endl;
    }
    std::cout << "\tChiSq = " << *(MinuitPolyFit.statistic) << std::endl;

    for (size_t i = 0; i < 3; ++i) {
        std::cout << "Polyfit params (no integral): " << MinuitPolyFitNoIntegral.fitParams.fitParams[i] << "+-"
                  << MinuitPolyFitNoIntegral.fitParams.fitParamErrors[i] << std::endl;
    }
    std::cout << "\tChiSq = " << *(MinuitPolyFitNoIntegral.statistic) << std::endl;

    // Plot fits to file
    CernFitter.plot->SetTitle("Compare Minuit and ROOT fitters;time/ns;DCS/CF ratio");
    CernFitter.plot->SetLineColor(kBlack);

    MinuitPolyFit.bestFitFunction->SetLineColor(kBlue);

    MinuitPolyFitNoIntegral.bestFitFunction->SetLineColor(6);
    MinuitPolyFitNoIntegral.bestFitFunction->SetLineStyle(9);
    MinuitPolyFitNoIntegral.bestFitFunction->SetLineWidth(3);

    PhysFitter.bestFitFunction->SetLineColor(kGreen);
    PhysFitter.bestFitFunction->SetLineStyle(kDashed);

    const util::LegendParams_t legendParams = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = "Compare fitters"};
    const std::vector<std::string> legendLabels{"Root best fit (red)",
                                                "Minuit polynomial best fit",
                                                "Minuit polynomial best fit (no integral)",
                                                "Minuit best fit"};

    util::saveObjectsToFile<TGraph>(std::vector<TObject*>{CernFitter.plot.get(),
                                                          MinuitPolyFit.bestFitFunction.get(),
                                                          MinuitPolyFitNoIntegral.bestFitFunction.get(),
                                                          PhysFitter.bestFitFunction.get()},
                                    std::vector<std::string>{"AP", "CSAME", "CSAME", "CSAME"},
                                    legendLabels,
                                    "compareMinuitRootPlots.pdf",
                                    legendParams);
}

int main()
{
    compareRootMinuit();
    return 0;
}
