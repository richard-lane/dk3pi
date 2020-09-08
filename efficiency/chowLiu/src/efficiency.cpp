#include <TH1D.h>
#include <TH2D.h>

#include "efficiency.h"
#include "efficiencyUtil.h"
#include "util.h"

ChowLiuEfficiency::ChowLiuEfficiency(const PhspBins& bins, const size_t root)
    : _detectedEvents(bins, "detectedEvents"), _generatedEvents(bins, "generatedEvents"), _root(root)
{
    ;
}

void ChowLiuEfficiency::addMCEvent(const PhspPoint& point)
{
    _detectedEvents.binPoint(point);
}

void ChowLiuEfficiency::addGeneratedEvent(const PhspPoint& point)
{
    _generatedEvents.binPoint(point);
}

void ChowLiuEfficiency::efficiencyParametrisation(void)
{
    _detectedEvents.makeApproximation();
    _generatedEvents.makeApproximation();

    _avgEfficiency     = (double)_detectedEvents.getNumPoints() / _generatedEvents.getNumPoints();
    _approximationMade = true;
}

double ChowLiuEfficiency::value(const PhspPoint& point) const
{
    if (!_approximationMade) {
        throw ChowLiu::ApproximationNotYetMade();
    }

    // The efficiency is p(detecting an event at x) / p(an event occurring at x)
    // Denominator comes directly from the generated-data histograms
    // Numerator comes from p(detected event being at x) * avg. efficiency
    double efficiency = _detectedEvents.value(point) * _avgEfficiency / _generatedEvents.value(point);

    // Just make a cursory check that things are sensible
    // if (efficiency > 1 || efficiency <= 0.0) {
    //    std::cout << _detectedEvents.value(point) << std::endl;
    //    std::cout << _generatedEvents.value(point) << std::endl;
    //    throw ChowLiu::BadProbability(efficiency, point);
    //}

    return efficiency;
}
