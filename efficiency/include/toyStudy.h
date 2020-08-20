#ifndef TOYSTUDY_H
#define TOYSTUDY_H

#include <functional>
#include <random>
#include <vector>

#include "efficiencyUtil.h"

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

#endif // TOYSTUDY_H
