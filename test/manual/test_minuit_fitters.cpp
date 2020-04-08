#include <iostream>
#include <vector>

#include "DecaySimulator.h"
#include "Fitter.h"
#include "PullStudyHelpers.h"
#include "RatioCalculator.h"
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

    // Create three fitters
    FitData_t  MyFitData         = FitData(MyRatios.binCentres, MyRatios.binWidths, MyRatios.ratio, MyRatios.error);
    RootFitter CernFitter        = RootFitter(MyFitData);
    Fitter     MinuitPolyFit     = Fitter(MyFitData);
    Fitter     MinuitDetailedFit = Fitter(MyFitData);

    // Perform fits
    std::vector<double> initialParameterGuess{0.02, 1.0, 100.0};
    std::vector<double> initialErrorsGuess{0.01, 1.0, 100.0};
    CernFitter.fit(0, maxTime, "Q");
    MinuitPolyFit.fitUsingMinuit(initialParameterGuess, initialErrorsGuess, ChiSquared);

    std::vector<double> initialParamGuess{phaseSpaceParams.x,
                                          phaseSpaceParams.y,
                                          phaseSpaceParams.r,
                                          phaseSpaceParams.z_im,
                                          phaseSpaceParams.z_re,
                                          phaseSpaceParams.width};
    std::vector<double> initialErrGuess{1, 1, 1, 1, 1, 1};
    MinuitDetailedFit.detailedFitUsingMinuit(initialParamGuess, initialErrGuess, ChiSquared);

    // Print fit parameters to console
    for (size_t i = 0; i < 3; ++i) {
        std::cout << "ROOT fit params: " << CernFitter.fitParams.fitParams[i] << "+-"
                  << CernFitter.fitParams.fitParamErrors[i] << std::endl;
    }
    std::cout << "\tChiSq = " << *(CernFitter.statistic) << std::endl;
    for (size_t i = 0; i < 3; ++i) {
        std::cout << "Chisq fit params: " << MinuitPolyFit.fitParams.fitParams[i] << "+-"
                  << MinuitPolyFit.fitParams.fitParamErrors[i] << std::endl;
    }
    std::cout << "\tChiSq = " << *(MinuitPolyFit.statistic) << std::endl;

    // Plot fits to file
    const util::LegendParams_t legendParams = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = "Compare fitters"};
    const std::vector<std::string> legendLabels{"Root best fit", "Minuit polynomial best fit", "Minuit best fit"};
    CernFitter.plot->SetTitle("Compare Minuit and ROOT fitters;time/ns;DCS/CF ratio");
    CernFitter.plot->SetLineColor(kBlack);
    MinuitPolyFit.bestFitPlot->SetLineColor(kBlue);
    MinuitDetailedFit.bestFitPlot->SetLineColor(kGreen);
    MinuitDetailedFit.bestFitPlot->SetLineStyle(kDashed);
    util::saveObjectsToFile<TGraph>(std::vector<TObject*>{CernFitter.plot.get(),
                                                          MinuitPolyFit.bestFitPlot.get(),
                                                          MinuitDetailedFit.bestFitPlot.get()},
                                    std::vector<std::string>{"AP", "CSAME", "CSAME"},
                                    legendLabels,
                                    "compareMinuitRootPlots.pdf",
                                    legendParams);
}

int main()
{
    compareRootMinuit();
    return 0;
}
