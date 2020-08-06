#include "efficiency.h"

Efficiency::Efficiency(const std::vector<dDecay_t>& detectedEvents,
                       const std::vector<dDecay_t>& generatedEvents,
                       const PhspBins<double>&      bins)
    : _detectedEvents(detectedEvents), _generatedEvents(generatedEvents), _bins(bins)
{
    ;
}

double Efficiency::value(const dDecay_t& event) const
{
    return event.dParams.energy;
}
