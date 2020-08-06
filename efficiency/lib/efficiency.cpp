#include "efficiency.h"
#include "efficiencyUtil.h"

Efficiency::Efficiency(const std::vector<dDecay_t>& detectedEvents,
                       const std::vector<dDecay_t>& generatedEvents,
                       const PhspBins<double>&      bins)
    : _detectedEvents(detectedEvents), _generatedEvents(generatedEvents), _bins(bins)
{
    // Find the invariant masses of our events
    _findInvariantMasses();
}

double Efficiency::value(const dDecay_t& event) const
{
    return event.dParams.energy;
    // Find the event's invariant mass

    // Find the m12, m23, m34 bin edges that are relevant for this event

    // Find how many generated events fall within this bin

    // Find how many detected events fall within this bin

    // return detected/generated
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
