#include "scriptUtils.h"

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

PhspPoint parametrisation(const dDecay_t& decay)
{
    // 5d phsp
    PhspPoint point = PhspPoint(5);

    // Use invariant masses m12, m23, m34, m123, m234
    point[0] = invariantMass({decay.kParams, decay.pi1Params});
    point[1] = invariantMass({decay.pi1Params, decay.pi2Params});
    point[2] = invariantMass({decay.pi2Params, decay.pi3Params});
    point[3] = invariantMass({decay.kParams, decay.pi1Params, decay.pi2Params});
    point[4] = invariantMass({decay.pi1Params, decay.pi2Params, decay.pi3Params});

    return point;
}

PhspBins findBins(void)
{
    std::cout << "Creating bins..." << std::flush;
    std::array<size_t, 5>                    numBins    = {100, 100, 100, 100, 100};
    std::array<std::pair<double, double>, 5> axisLimits = {std::make_pair(0.4, 1.6),
                                                           std::make_pair(0.2, 1.4),
                                                           std::make_pair(0.2, 1.4),
                                                           std::make_pair(0.8, 1.8),
                                                           std::make_pair(0.4, 1.8)};

    PhspBins Bins(5);
    for (size_t i = 0; i < Bins.size(); ++i) {
        Bins[i] = std::vector<double>(numBins[i] + 1);
        for (size_t j = 0; j <= numBins[i]; ++j) {
            Bins[i][j] = axisLimits[i].first + (axisLimits[i].second - axisLimits[i].first) * j / (numBins[i]);
        }
    }

    std::cout << "done" << std::endl;
    return Bins;
}