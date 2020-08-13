#include <TH1D.h>
#include <TH2D.h>

#include "efficiency.h"
#include "efficiencyUtil.h"
#include "util.h"

EfficiencyBinning::EfficiencyBinning(const PhspBins& bins) : _bins(bins)
{
    // Need to set the numbers of bins first
    for (size_t i = 0; i < _bins.size(); ++i) {
        _numBins[i] = _bins[i].size() - 1;
    }

    for (size_t i = 0; i < _bins.size(); ++i) {
        std::string title1d = "1dhist" + std::to_string(i);
        _1dhistograms[i]    = std::make_unique<TH1D>(title1d.c_str(), title1d.c_str(), _numBins[i], _bins[i].data());

        for (size_t j = _bins.size() - 1; j > i; --j) {
            std::string title2d                   = "2dhist(" + std::to_string(i) + "," + std::to_string(j) + ")";
            _2dhistograms[_indexConversion(i, j)] = std::make_unique<TH2D>(
                title2d.c_str(), title2d.c_str(), _numBins[i], _bins[i].data(), _numBins[j], _bins[j].data());
        }
    }
}

EfficiencyBinning::EfficiencyBinning(const PhspBins& bins, const PhspPoint& point) : EfficiencyBinning(bins)
{
    binPoint(point);
}

void EfficiencyBinning::binPoint(const PhspPoint& point)
{
    for (size_t i = 0; i < _bins.size(); ++i) {
        _1dhistograms[i]->Fill(point[i]);

        for (size_t j = _bins.size() - 1; j > i; --j) {
            _2dhistograms[_indexConversion(i, j)]->Fill(point[i], point[j]);
        }
    }
}

const TH1D EfficiencyBinning::get1dhistogram(const unsigned short i) const
{
    if (i > 4) {
        throw HistogramNotFound();
    }
    return *_1dhistograms[i];
}

const TH2D EfficiencyBinning::get2dhistogram(const unsigned short i, const unsigned short j) const
{
    return *_2dhistograms[_indexConversion(i, j)];
}

unsigned short EfficiencyBinning::_indexConversion(const unsigned short i, const unsigned short j) const
{
    if (i == j) {
        std::cerr << "No 2d histogram of the a variable against itself exists" << std::endl;
        throw HistogramNotFound();
    }

    // convert i, j -> smaller, larger
    unsigned short smaller = i < j ? i : j;
    unsigned short larger  = i > j ? i : j;

    if (larger > 4) {
        throw HistogramNotFound();
    }

    // We store our 2d histograms in a 1d array
    // The array looks like {i,j} = {04, 03, 02, 01, 14, 13, 12, 11, 24, 23, 34}
    return (9 * smaller - smaller * smaller) / 2 + 4 - larger;
}

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

    size_t numEvents = xHist->GetEntries(); // For some reason histogram2d->GetEntries() doesn't give the right
                                            // number of events so just do this
    for (size_t xBin = 1; xBin <= numBinsX; ++xBin) {
        for (size_t yBin = 1; yBin <= numBinsY; ++yBin) {
            size_t binContent = histogram2d->GetBinContent(xBin, yBin);
            if (binContent == 0) {
                // I don't like breaking the control flow here but i'm tired and can't immediately think of
                // something better. maybe an if/else
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
