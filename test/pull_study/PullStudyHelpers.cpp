#ifndef PULL_STUDY_HELPERS_CPP
#define PULL_STUDY_HELPERS_CPP

#include <iostream>
#include <string>
#include <vector>

#include "TH1D.h"

#include "DecaySimulator.h"
#include "physics.h"
#include "util.h"

namespace PullStudyHelpers
{

std::pair<double, double> meanAndStdDev(const std::vector<double>& v)
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

    TH1D* MyGraph = new TH1D(title.c_str(), title.c_str(), 200, xMin, xMax);

    MyGraph->FillN(nExperiments, parameter.data(), 0);
    MyGraph->SetTitle((title + ";Normalised Pull;Count").c_str());

    util::saveObjectToFile(MyGraph, (title + ".pdf").c_str());

    std::cout << title + " mean:\t\t" + MyGraph->GetMean() << std::endl;
    std::cout << title + " std dev:\t" + MyGraph->GetStdDev() << std::endl;
    delete MyGraph;
}

/*
 * Plot a histogram from a vector
 */
void plotHist(const std::vector<double>& vector, const size_t numBins, const std::string& name)
{
    assert((numBins != 0));
    double min = *(std::min_element(vector.begin(), vector.end()));
    double max = *(std::max_element(vector.begin(), vector.end()));

    double              binWidth = (max - min) / numBins;
    std::vector<double> binLimits(numBins + 1, -1);

    binLimits[0]       = min * 0.99;
    binLimits[numBins] = max * 1.01;
    for (size_t i = 1; i < numBins + 1; i++) {
        binLimits[i] = binWidth * i + min;
    }

    TH1D* hist = new TH1D(name.c_str(), name.c_str(), numBins, binLimits.data());
    hist->FillN(vector.size(), vector.data(), 0);

    util::saveObjectToFile(hist, name + ".pdf");
    delete hist;
}

double numDCSDecays(const size_t numCFDecays, const DecayParams_t& phaseSpaceParams, double maxTime)
{
    // Our formula is numDcs = numCf * (DCS integral / CF integral), where we integrate over all allowed times
    double dcsIntegral =
        Phys::dcsIntegralWithEfficiency(0, maxTime, util::expectedParams(phaseSpaceParams), phaseSpaceParams.width);
    double cfIntegral = Phys::cfIntegralWithEfficiency(0, maxTime, phaseSpaceParams.width);

    return numCFDecays * dcsIntegral / cfIntegral;
}

} // namespace PullStudyHelpers

#endif // PULL_STUDY_HELPERS_CPP
