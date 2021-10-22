#include <assert.h>
#include <stdlib.h>
#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "CleoCombinationFitter.h"
#include "bindings.h"

std::vector<double> cleoZScan(const std::vector<double>&  reZVals,
                             const std::vector<double>&   imZVals,
                             const std::vector<double>&   rsDecayTimes,
                             const std::vector<double>&   rsWeights,
                             const std::vector<double>&   wsDecayTimes,
                             const std::vector<double>&   wsWeights,
                             const std::vector<double>&   binLimits,
                             const std::array<double, 6>& initialVals,
                             const std::array<double, 6>& initialErrs,
                             const int                    binNumber)
{
    auto phspBin = static_cast<CLEO::Bin>(binNumber);

    const auto nReVals = reZVals.size();
    const auto nImVals = imZVals.size();

    const unsigned nVals = nReVals * nImVals;
    std::vector<double> likelihoods(nVals);

    // Create a fitter
    CharmFitter::CLEOCombinationFitter Fitter(binLimits, initialVals, initialErrs, phspBin);
    Fitter.fixParameters(std::array<std::string, 3>{"width", "z_im", "z_re"});

    Fitter.addRSPoints(rsDecayTimes, rsWeights);
    Fitter.addWSPoints(wsDecayTimes, wsWeights);

    // Efficiency = 1 since our weights account for the decay time efficiency...
    const auto efficiency = []([[maybe_unused]] const double x) { return 1; };

    for (size_t rIndex=0; rIndex < nReVals; ++rIndex) {
        for (size_t iIndex=0; iIndex < nImVals; ++iIndex) {
            Fitter.setParameter("z_re", reZVals[rIndex]);
            Fitter.setParameter("z_im", imZVals[iIndex]);

            const unsigned index = iIndex + rIndex * nImVals;
            likelihoods[index]   = Fitter.fit(efficiency).fitStatistic;
        }
    }
    return likelihoods;
}

