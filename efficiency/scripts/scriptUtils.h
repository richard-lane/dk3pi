#ifndef EFFICIENCY_SCRIPT_UTILS_H
#define EFFICIENCY_SCRIPT_UTILS_H

#include <iostream>
#include <memory>
#include <random>
#include <string>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLegend.h>

#include "ReadRoot.h"
#include "efficiency.h"
#include "efficiencyUtil.h"
#include "util.h"

struct EventDetectionProbNotNormalised : public std::exception {
    EventDetectionProbNotNormalised(const double prob)
        : _msg("Probability " + std::to_string(prob) + " not between 0 and 1")
    {
        ;
    }

    const char* what() const throw() { return _msg.c_str(); }

  private:
    const std::string _msg;
};

/*
 * Provide a random number generator, a vector of D decay events and a function that provides an event's detection
 * probability
 *
 * The probability of detecting an event should be a number between 0 and 1
 *
 * Removes the undetected events from the vector in-place.
 */
void applyEfficiency(std::mt19937* const                    generator,
                     const std::function<double(dDecay_t)>& eventDetectionProb,
                     std::vector<dDecay_t>&                 events);

/*
 * Return the 5d parametrisation of a phase space point
 */
PhspPoint parametrisation(const dDecay_t& decay);

inline double simpleEfficiency([[maybe_unused]] const dDecay_t& event)
{
    return 0.5;
}

/*
 * A nice efficiency that gets recovered well
 */
inline double niceEfficiency(const dDecay_t& event)
{
    return invariantMass({event.kParams, event.pi1Params}) / 3;
}

/*
 * A less-nice efficiency that doesn't get recovered as nicely
 */
inline double awkwardEfficiency(const dDecay_t& event)
{
    std::vector<double> params = parametrisation(event);
    double              e{1};
    for (int i = 0; i < 3; ++i) {
        e *= params[i] / 2;
    }
    return 5 * e;
}

/*
 * An efficiency on the total pT of the k
 */
inline double pTEfficiency(const dDecay_t& event)
{
    return std::sqrt(pT(event.kParams));
}

/*
 * Return some sensible looking bin limits
 */
PhspBins findBins(void);

#endif // EFFICIENCY_SCRIPT_UTILS_H
