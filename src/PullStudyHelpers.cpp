#ifndef PULL_STUDY_HELPERS_CPP
#define PULL_STUDY_HELPERS_CPP

#include <iostream>
#include <string>
#include <vector>

#include "TH1D.h"

#include "../lib/DecaySimulator.h"
#include "../lib/util.h"

namespace PullStudyHelpers
{
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

    util::saveToFile(MyGraph, (title + ".pdf").c_str());

    std::cout << title + " mean:\t\t" + MyGraph->GetMean() << std::endl;
    std::cout << title + " std dev:\t" + MyGraph->GetStdDev() << std::endl;
    delete MyGraph;
}

std::vector<double> expectedParams(const DecayParams_t &phaseSpaceParams)
{
    double expected_a = phaseSpaceParams.r * phaseSpaceParams.r;
    double expected_b = phaseSpaceParams.r *
                        (phaseSpaceParams.y * phaseSpaceParams.z_re + phaseSpaceParams.x * phaseSpaceParams.z_im) *
                        phaseSpaceParams.width;
    double expected_c = 0.25 * (std::pow(phaseSpaceParams.x, 2) + std::pow(phaseSpaceParams.y, 2)) *
                        std::pow(phaseSpaceParams.width, 2);

    return std::vector<double>{expected_a, expected_b, expected_c};
}

size_t numDCSDecays(const size_t numCFDecays, const DecayParams_t &phaseSpaceParams, double maxTime)
{
    // Our formula is prefactor * integral * numCFDecays
    double              exp        = std::exp(-1 * phaseSpaceParams.width * maxTime);
    std::vector<double> fit_params = expectedParams(phaseSpaceParams);
    double              a          = fit_params[0];
    double              b          = fit_params[1];
    double              c          = fit_params[2];

    double prefactor = phaseSpaceParams.width / (1 - exp * maxTime);

    // crikey
    // Put the integral (a + bt + ct^3)e^-width*t into WolframAlpha, and hopefully this is what you get
    // TODO UT PLEASE
    double integral =
        exp * ((-1 * phaseSpaceParams.width * (a * phaseSpaceParams.width + b * phaseSpaceParams.width * maxTime + b) -
                c * (std::pow(phaseSpaceParams.width * maxTime, 2) + 2 * phaseSpaceParams.width * maxTime + 2)) /
               std::pow(phaseSpaceParams.width, 3)) +
        (phaseSpaceParams.width * (a * phaseSpaceParams.width + b) + 2 * c) / std::pow(phaseSpaceParams.width, 3);

    return prefactor * integral * numCFDecays;
}

} // namespace PullStudyHelpers

#endif // PULL_STUDY_HELPERS_CPP
