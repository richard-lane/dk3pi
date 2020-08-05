#include <random>

#include "toyStudy.h"

void applyEfficiency(std::mt19937* const                    generator,
                     const std::function<double(dDecay_t)>& eventDetectionProb,
                     std::vector<dDecay_t>&                 events)
{
    std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);

    // Lambda that we'll use to decide whether to remove our event
    auto removeEvent = [&](const dDecay_t& event) {
        double detectionProb = eventDetectionProb(event);
        if (detectionProb < 0 || detectionProb > 1) {
            throw EventDetectionProbNotNormalised(detectionProb);
        }
        return uniformDistribution(*generator) > detectionProb;
    };

    // Move the items to remove to the end of the vector and erase them
    auto it = std::remove_if(events.begin(), events.end(), removeEvent);
    events.erase(it, events.end());
}
