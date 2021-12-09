#include <assert.h>
#include <stdlib.h>
#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "CleoCombinationFitter.h"
#include "bindings.h"

std::vector<double> cleoLikelihoods(const std::vector<double>&   reZVals,
                                    const std::vector<double>&   imZVals,
                                    const std::array<double, 6>& decayParams,
                                    const int                    binNumber)
{
    auto phspBin = static_cast<CLEO::Bin>(binNumber);

    const auto nReVals = reZVals.size();
    const auto nImVals = imZVals.size();

    const unsigned nVals = nReVals * nImVals;
    std::vector<double> likelihoods(nVals);

    auto params = FitterUtil::DecayParams_t{decayParams[0],
                                       decayParams[1],
                                       decayParams[2],
                                       decayParams[3],
                                       decayParams[4],
                                       decayParams[5]};

    for (size_t i=0; i<nReVals; ++i) {
        for (size_t j=0; j<nImVals; ++j) {
            // Build decay param struct
            params.z_re = reZVals[i];
            params.z_im = imZVals[j];

            const unsigned index = j + i * nImVals;  // TODO this might be wrong
            likelihoods[index] = CLEO::cleoLikelihood(phspBin, params);
        }
    }

    return likelihoods;
}

static std::vector<double> _do_scan(CharmFitter::CharmFitterBase* Fitter,
                                    const std::vector<double>& reZVals,
                                    const std::vector<double>& imZVals) {
    // Efficiency = 1 since our weights account for the decay time efficiency...
    const auto efficiency = []([[maybe_unused]] const double x) { return 1.0; };

    const auto nReVals = reZVals.size();
    const auto nImVals = imZVals.size();
    const unsigned nVals = nReVals * nImVals;

    std::vector<double> likelihoods(nVals, 0);
    for (size_t rIndex=0; rIndex < nReVals; ++rIndex) {
        for (size_t iIndex=0; iIndex < nImVals; ++iIndex) {
            Fitter->setParameter("z_re", reZVals[rIndex]);
            Fitter->setParameter("z_im", imZVals[iIndex]);

            auto     res   = Fitter->fit(efficiency);
            likelihoods[iIndex + rIndex * nImVals]   = res.fitStatistic;
        }
    }
    return likelihoods;
}

std::vector<double> charmFitLikelihoods(const std::vector<double>&   reZVals,
                                        const std::vector<double>&   imZVals,
                                        const std::vector<double>&   rsDecayTimes,
                                        const std::vector<double>&   rsWeights,
                                        const std::vector<double>&   wsDecayTimes,
                                        const std::vector<double>&   wsWeights,
                                        const std::vector<double>&   binLimits,
                                        const std::array<double, 6>& initialVals,
                                        const std::array<double, 6>& initialErrs) {
    // Create a fitter
    auto Fitter = CharmFitter::ConstrainedFitter(binLimits, initialVals, initialErrs);
    Fitter.fixParameters(std::array<std::string, 3>{"width", "z_im", "z_re"});

    Fitter.addRSPoints(rsDecayTimes, rsWeights);
    Fitter.addWSPoints(wsDecayTimes, wsWeights);

    return _do_scan(&Fitter, reZVals, imZVals);
}


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

    // Create a fitter
    auto Fitter = CharmFitter::CLEOCombinationFitter(binLimits, initialVals, initialErrs, phspBin);
    Fitter.fixParameters(std::array<std::string, 3>{"width", "z_im", "z_re"});

    Fitter.addRSPoints(rsDecayTimes, rsWeights);
    Fitter.addWSPoints(wsDecayTimes, wsWeights);

    auto likelihoods = _do_scan(&Fitter, reZVals, imZVals);

    // Find the maximum value inside the allowed region
    double max = 0;
    for (size_t i=0; i < reZVals.size(); ++i) {
        for (size_t j=0; j < imZVals.size(); ++j) {
            if ((reZVals[i] * reZVals[i] + imZVals[j] * imZVals[j]) < 1.0 && likelihoods[j + i * imZVals.size()] > max) {
                max = likelihoods[j + i * imZVals.size()];
            }
        }
    }

    // Set un allowed values to this
    for (size_t i=0; i < reZVals.size(); ++i) {
        for (size_t j=0; j < imZVals.size(); ++j) {
            if ((reZVals[i] * reZVals[i] + imZVals[j] * imZVals[j]) > 1.0) {
                likelihoods[j + i * imZVals.size()] = max;
            }
        }
    }

    return likelihoods;
}

