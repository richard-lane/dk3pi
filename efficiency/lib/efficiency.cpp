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

double Efficiency::value(const std::vector<double>& invMasses) const
{

    // Find the m12, m23, m34, m123, m234 bin edges that are relevant for this event
    std::vector<std::pair<double, double>> eventBinEdges(5);
    for (size_t i = 0; i < invMasses.size(); ++i) {
        // Find what bin the event belongs in by binning a single event and finding where it ends up
        std::vector<size_t> binLocation    = util::binVector(std::vector<double>{invMasses}, _bins[i]);
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

double entropy(const TH1D* const hist)
{
    size_t totalEntries = hist->GetEntries();
    double ent{0.0}; // entropy

    for (size_t bin = 1; bin <= (size_t)hist->GetNbinsX(); ++bin) { // ROOT bin indexing starts at 1
        double prob = (double)hist->GetBinContent(bin) / totalEntries;
        ent -= prob * std::log(prob);
    }
    return ent;
}

double mutual_info(const TH2D* const histogram2d)
{
    double info{0.0};

    // A better implementation probably wouldn't construct an entire TH1 just to find the sum along an axis
    TH1D* xHist = histogram2d->ProjectionX();
    TH1D* yHist = histogram2d->ProjectionY();

    size_t numBinsX = xHist->GetNbinsX();
    size_t numBinsY = yHist->GetNbinsX();

    size_t numEvents = xHist->GetEntries(); // For some reason histogram2d->GetEntries() doesn't give the right number
                                            // of events so just do this
    for (size_t xBin = 1; xBin <= numBinsX; ++xBin) {
        for (size_t yBin = 1; yBin <= numBinsY; ++yBin) {
            size_t binContent = histogram2d->GetBinContent(xBin, yBin);
            if (binContent == 0) {
                // I don't like breaking the control flow here but i'm tired and can't immediately think of something
                // better. maybe an if/else
                continue;
            }

            size_t xBinContent = xHist->GetBinContent(xBin);
            size_t yBinContent = yHist->GetBinContent(yBin);
            // These casts aren't all necessary
            double thingInBrackets =
                (double)numEvents * (double)binContent / ((double)xBinContent * (double)yBinContent);
            info += ((double)binContent / (double)numEvents) * std::abs(std::log(thingInBrackets));
        }
    }

    // Normalise our mutual information
    double normalisation = 2 / (entropy(xHist) + entropy(yHist));

    return normalisation * info;
}
