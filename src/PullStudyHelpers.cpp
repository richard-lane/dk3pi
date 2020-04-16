#ifndef PULL_STUDY_HELPERS_CPP
#define PULL_STUDY_HELPERS_CPP

#include <boost/math/quadrature/trapezoidal.hpp>
#include <iostream>
#include <string>
#include <vector>

#include "TH1D.h"

#include "../lib/DecaySimulator.h"
#include "../lib/util.h"

namespace PullStudyHelpers
{

std::vector<double> exponentialBinLimits(const double maxTime, const double decayConstant, const size_t numBins)
{
    std::vector<double> binLimits{};
    for (size_t i = 0; i <= numBins; ++i) {
        double x = (double)i / numBins;
        double z = 1 - std::exp(-1 * decayConstant * maxTime);
        binLimits.push_back((-1 / decayConstant) * std::log(1 - z * x));
    }

    return binLimits;
}

std::pair<double, double> meanAndStdDev(const std::vector<double> &v)
{

    double sum  = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();

    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(), std::bind2nd(std::minus<double>(), mean));
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev  = std::sqrt(sq_sum / v.size());

    return std::make_pair(mean, stdev);
}

void plot_parameter_distribution(std::string         title,
                                 std::vector<double> parameter,
                                 size_t              nExperiments,
                                 double              expectedMean,
                                 double              expectedSigma)
{
    // Define axis limits
    double xMin = expectedMean - 5 * expectedSigma;
    double xMax = expectedMean - 5 * expectedSigma;

    TH1D *MyGraph = new TH1D(title.c_str(), title.c_str(), 200, xMin, xMax);

    MyGraph->FillN(nExperiments, parameter.data(), 0);
    MyGraph->SetTitle((title + ";Normalised Pull;Count").c_str());

    util::saveObjectToFile(MyGraph, (title + ".pdf").c_str());

    std::cout << title + " mean:\t\t" + MyGraph->GetMean() << std::endl;
    std::cout << title + " std dev:\t" + MyGraph->GetStdDev() << std::endl;
    delete MyGraph;
}

double numDCSDecays(const size_t numCFDecays, const DecayParams_t &phaseSpaceParams, double maxTime)
{
    // Our formula is prefactor * integral * numCFDecays
    double              exp        = std::exp(-1 * phaseSpaceParams.width * maxTime);
    std::vector<double> fit_params = util::expectedParams(phaseSpaceParams);
    double              a          = fit_params[0];
    double              b          = fit_params[1];
    double              c          = fit_params[2];

    double prefactor = phaseSpaceParams.width / (1 - exp);

    auto   f        = [&](double x) { return (a + b * x + c * x * x) * std::exp(-phaseSpaceParams.width * x); };
    double integral = boost::math::quadrature::trapezoidal(f, 0.0, maxTime, 1e-10, 20);
    return prefactor * integral * numCFDecays;
}

} // namespace PullStudyHelpers

#endif // PULL_STUDY_HELPERS_CPP
