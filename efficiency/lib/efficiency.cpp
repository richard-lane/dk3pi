#include "efficiency.h"
#include "efficiencyUtil.h"
#include "util.h"

Efficiency::Efficiency(const std::vector<dDecay_t>&            detectedEvents,
                       const std::vector<dDecay_t>&            generatedEvents,
                       const std::vector<std::vector<double>>& bins)
    : _detectedEvents(detectedEvents), _generatedEvents(generatedEvents), _bins(bins)
{
    // Find the invariant masses of our events
    _findInvariantMasses();
}

double Efficiency::value(const dDecay_t& event) const
{
    return event.dParams.energy;
    // Find the event's invariant masses
    std::vector<double> invariantMasses = event2invariantMasses(event);

    // Find the m12, m23, m34, m123, m234 bin edges that are relevant for this event
    std::vector<std::pair<double, double>> eventBinEdges(5);
    for (size_t i = 0; i < invariantMasses.size(); ++i) {
        // Find what bin the event belongs in by binning a single event and finding where it ends up
        std::vector<size_t> binLocation    = util::binVector(std::vector<double>{invariantMasses}, _bins[i]);
        auto                lowBinIterator = std::find(binLocation.begin(), binLocation.end(), 1);

        eventBinEdges[i].first  = *lowBinIterator;
        eventBinEdges[i].second = *(lowBinIterator + 1);
    }

    // Find how many generated events fall within this bin
    size_t numGenEventsInBin{0};
    for (auto genEvent : _generatedEvents) {
        std::vector<double> genInvMasses = event2invariantMasses(genEvent);
        for (size_t i = 0; i < genInvMasses.size(); ++i) {
            if (genInvMasses[i] < eventBinEdges[i].first || genInvMasses[i] > eventBinEdges[i].second) {
                break; // I don't like this
            }
            numGenEventsInBin++;
        }
    }

    // If 0 generated events, return 0
    if (numGenEventsInBin == 0) {
        return 0;
    }

    // Find how many detected events fall within this bin
    size_t numDetectedEventsInBin{0};
    for (auto detectedEvent : _detectedEvents) {
        std::vector<double> detectedInvMasses = event2invariantMasses(detectedEvent);
        for (size_t i = 0; i < detectedInvMasses.size(); ++i) {
            if (detectedInvMasses[i] < eventBinEdges[i].first || detectedInvMasses[i] > eventBinEdges[i].second) {
                break; // I don't like this
            }
            numDetectedEventsInBin++;
        }
    }

    // return detected/generated
    return (double)numDetectedEventsInBin / (double)numGenEventsInBin;
}

void Efficiency::_findInvariantMasses(void)
{
    _detectedInvMasses  = std::vector<std::vector<double>>(_detectedEvents.size());
    _generatedInvMasses = std::vector<std::vector<double>>(_generatedEvents.size());

    for (size_t i = 0; i < _detectedEvents.size(); ++i) {
        _detectedInvMasses[i] = event2invariantMasses(_detectedEvents[i]);
    }

    for (size_t i = 0; i < _generatedEvents.size(); ++i) {
        _generatedInvMasses[i] = event2invariantMasses(_generatedEvents[i]);
    }
}
