/*
 * Input ideal, noiseless ratios to our fitter and hope that it outputs the right fit parameters
 */

#include <iostream>
#include <random>
#include <vector>

#include "DecaySimulator.h"
#include "Fitter.h"
#include "PullStudyHelpers.h"

#include "TGraph.h"

/*
 * Find the expected ratio at a given time
 */
inline double ratio(const double a, const double b, const double c, const double time)
{
    return a + b * time + c * time * time;
}

/*
 * Create an idealised set of ratios and times with the provided errors, fit them with the Fitter() class and return a
 * vector of fit parameters.
 *
 * Takes loads of parameters because i cba to make it better
 *
 * returns the fitter cus why not
 */
Fitter performFit(size_t numTimeBins, double maxTime, double a, double b, double c, double error)
{
    std::random_device               randomEngine;
    std::mt19937                     mt(randomEngine());
    std::normal_distribution<double> distribution(0.0, 1.0);

    double              timeBinWidth = maxTime / numTimeBins;
    std::vector<double> times        = std::vector<double>(numTimeBins, -1);
    std::vector<double> ratios       = std::vector<double>(numTimeBins, -1);
    std::vector<double> ratioErrors  = std::vector<double>(numTimeBins, -1);

    // Create idealised plot
    for (size_t i = 0; i < numTimeBins; ++i) {
        double time      = i * timeBinWidth;
        double thisRatio = ratio(a, b, c, time);
        times[i]         = time;
        ratios[i]        = thisRatio * (1 + error * distribution(mt));
        ratioErrors[i]   = thisRatio * error;
    }

    // Fit our idealised plot
    FitData_t MyFitData = FitData(times, std::vector<double>(numTimeBins, timeBinWidth), ratios, ratioErrors);
    Fitter    MyFitter(MyFitData);
    MyFitter.fitUsingMinuit(std::vector<double>{a, b, c}, std::vector<double>{1, 1, 1}, ChiSquared);
    // MyFitter.saveFitPlot("Generated DCS/CF ratios", "baz.pdf");

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

    double maxTime     = 0.005;
    size_t numTimeBins = 20;

    // We want to add random gaussian noise to our ratios
    // The error is calculated from error * random number
    double error = 1;

    size_t              numExperiments = 100000;
    std::vector<double> aPull          = std::vector<double>(numExperiments, std::nan("-1"));
    std::vector<double> bPull          = std::vector<double>(numExperiments, std::nan("-1"));
    std::vector<double> cPull          = std::vector<double>(numExperiments, std::nan("-1"));

    for (size_t i = 0; i < numExperiments; ++i) {
        Fitter              MyFitter  = performFit(numTimeBins, maxTime, a, b, c, error);
        std::vector<double> outParams = MyFitter.fitParams.fitParams;
        std::vector<double> outErrors = MyFitter.fitParams.fitParamErrors;
        aPull[i]                      = (outParams[0] - a) / outErrors[0];
        bPull[i]                      = (outParams[1] - b) / outErrors[1];
        cPull[i]                      = (outParams[2] - c) / outErrors[2];
    }

    PullStudyHelpers::plot_parameter_distribution("a", aPull, numExperiments);
    PullStudyHelpers::plot_parameter_distribution("b", bPull, numExperiments);
    PullStudyHelpers::plot_parameter_distribution("c", cPull, numExperiments);

    return 0;
}
