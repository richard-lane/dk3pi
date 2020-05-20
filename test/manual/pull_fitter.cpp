/*
 * Input ideal, noiseless ratios to our fitter and hope that it outputs the right fit parameters
 */

#include <iostream>
#include <random>
#include <vector>

#include "../pull_study/DecaySimulator.h"
#include "../pull_study/PullStudyHelpers.h"
#include "fitter/MinuitPolynomialFitter.h"

#include "TGraph.h"

#include "testutil.h"

#include <boost/progress.hpp>

/*
 * Create an idealised set of ratios and times with the provided errors, fit them with the Fitter() class and return a
 * vector of fit parameters.
 *
 * Takes loads of parameters because i cba to make it better
 *
 * returns the fitter cus why not
 */
MinuitPolynomialFitter
performFit(size_t numTimeBins, double maxTime, double a, double b, double c, double error, const double width)
{
    std::random_device               randomEngine;
    std::mt19937                     mt(randomEngine());
    std::normal_distribution<double> distribution(0.0, 1.0);

    std::vector<double> timeBinWidths = std::vector<double>(numTimeBins, -1);
    std::vector<double> times         = std::vector<double>(numTimeBins, -1);
    std::vector<double> ratioErrors   = std::vector<double>(numTimeBins, -1);

    // Create exponential bins
    std::vector<double> timeBinLimits = util::exponentialBinLimits(maxTime, width, numTimeBins);

    // Create idealised plot
    for (size_t i = 0; i < numTimeBins; ++i) {
        times[i]         = (timeBinLimits[i] + timeBinLimits[i + 1]) / 2;
        timeBinWidths[i] = (timeBinLimits[i + 1] - timeBinLimits[i]);
    }

    std::vector<double> ratios = idealRatios(times, error, a, b, c);
    for (size_t i = 0; i < numTimeBins; ++i) {
        ratioErrors[i] = ratio(a, b, c, times[i]) * error;
    }

    // Fit our idealised plot
    FitData_t              MyFitData = FitData(times, timeBinWidths, ratios, ratioErrors);
    MinuitPolynomialFitter MyFitter(MyFitData);
    MyFitter.setPolynomialParams(std::vector<double>{a, b, c}, std::vector<double>{1, 1, 1});
    MyFitter.fit();
    // const util::LegendParams_t legendParams = {
    //     .x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = "Compare fitters "};
    // MyFitter.saveFitPlot(" Generated DCS / CF ratios ", " baz.pdf ", &legendParams);

    return MyFitter;
}

int main()
{
    // Create an idealised plot of our ratios
    // Find the parameters we expect
    DecayParams_t DecayParams = {
        .x     = 0.004,
        .y     = 0.007,
        .r     = 0.05,
        .z_im  = -0.3,
        .z_re  = 0.8,
        .width = 2500.0,
    };

    std::vector<double> params = util::expectedParams(DecayParams);
    double              a      = params[0];
    double              b      = params[1];
    double              c      = params[2];

    double maxTime     = 0.002;
    size_t numTimeBins = 50;

    // We want to add random gaussian noise to our ratios
    // The error is calculated from error * random number
    double error = 1;

    size_t              numExperiments = 10000;
    std::vector<double> aPull          = std::vector<double>(numExperiments, std::nan("-1"));
    std::vector<double> bPull          = std::vector<double>(numExperiments, std::nan("-1"));
    std::vector<double> cPull          = std::vector<double>(numExperiments, std::nan("-1"));

    boost::progress_display showProgress(numExperiments);

    for (size_t i = 0; i < numExperiments; ++i) {
        MinuitPolynomialFitter MyFitter  = performFit(numTimeBins, maxTime, a, b, c, error, DecayParams.width);
        std::vector<double>    outParams = MyFitter.fitParams.fitParams;
        std::vector<double>    outErrors = MyFitter.fitParams.fitParamErrors;
        aPull[i]                         = (outParams[0] - a) / outErrors[0];
        bPull[i]                         = (outParams[1] - b) / outErrors[1];
        cPull[i]                         = (outParams[2] - c) / outErrors[2];
        ++showProgress;
    }

    PullStudyHelpers::plot_parameter_distribution("a", aPull, numExperiments);
    PullStudyHelpers::plot_parameter_distribution("b", bPull, numExperiments);
    PullStudyHelpers::plot_parameter_distribution("c", cPull, numExperiments);

    return 0;
}
