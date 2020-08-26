#include <TH1D.h>
#include <TH2D.h>

#include "efficiency.h"
#include "efficiencyUtil.h"
#include "util.h"

ChowLiuEfficiency::ChowLiuEfficiency(const PhspBins& bins, const size_t root)
    : _bins(bins), _root(root), _dimensionality(_bins.size())
{
    // Initialise the classes used to hold our data
    _detectedEvents = std::make_unique<ChowLiu::HistogramProjections>(ChowLiu::HistogramProjections(_bins, "detected"));
    _generatedEvents =
        std::make_unique<ChowLiu::HistogramProjections>(ChowLiu::HistogramProjections(_bins, "generated"));

    // Initialise the graph that we'll use to work out the best approximation
    _graph = std::make_unique<Graph>(_dimensionality);
}

void ChowLiuEfficiency::addMCEvent(const PhspPoint& point)
{
    _detectedEvents->binPoint(point);
}

void ChowLiuEfficiency::addGeneratedEvent(const PhspPoint& point)
{
    _generatedEvents->binPoint(point);
}

void ChowLiuEfficiency::efficiencyParametrisation(void)
{
    // Find the histograms for detected/generated events
    _ratio = ChowLiu::HistogramProjections(*_detectedEvents / *_generatedEvents);

    // Find the mutual information for each edge, add it to the graph
    for (size_t outNode = 0; outNode < _dimensionality; ++outNode) {
        for (size_t inNode = outNode + 1; inNode < _dimensionality; ++inNode) {
            TH2D   hist   = _ratio.get2dhistogram(inNode, outNode);
            double weight = ChowLiu::mutual_info(&hist);
            _graph->addEdge(outNode, inNode, weight);
        }
    }

    // Find the directed MST and store it
    _directedTree = inTree(_root, _graph->getMaxSpanningTree());

    // Find what variables our histograms are in
    for (auto edges : _directedTree) {
        for (Edge edge : edges) {
            _2dHistVars.push_back(std::make_pair(edge.from(), edge.to()));
        }
    }

    // Populate our histograms
    for (size_t i = 0; i < _dimensionality; ++i) {
        _hists1d.push_back(std::make_unique<TH1D>(_ratio.get1dhistogram(i)));
    }

    for (auto pair : _2dHistVars) {
        // The HistogramProjections class only stores 2d histograms (i, j) for i<j
        // We may need to swap the axes of our histogram if we want (j, i) for our parametrisation
        TH2D hist = pair.first > pair.second ? ChowLiu::swapAxes(_ratio.get2dhistogram(pair.first, pair.second))
                                             : _ratio.get2dhistogram(pair.first, pair.second);
        _hists2d.push_back(std::make_unique<TH2D>(hist));
    }

    _approximationMade = true;
}

double ChowLiuEfficiency::value(const PhspPoint& point) const
{
    if (!_approximationMade) {
        throw ChowLiu::ApproximationNotYetMade();
    }
    // Set the value to the 1d prob
    double prob1d = _hists1d[_root]->GetBinContent(_hists1d[_root]->FindBin(point[_root]));
    double prob{prob1d};

    // Iterate over 2d histograms (there are d-1 of them), finding the conditional probabilities
    for (size_t i = 0; i < _dimensionality - 1; ++i) {
        // The 2d histograms we stored are (x, y) for p(x|y)
        // Conditional efficiency is e(x & y)/e(y)

        // The x and y bins our point is in on our 2d histogram
        size_t xBin = _hists2d[i]->GetXaxis()->FindBin(point[_2dHistVars[i].first]);
        size_t yBin = _hists2d[i]->GetYaxis()->FindBin(point[_2dHistVars[i].second]);

        // Event detection prob for x&y
        double pXAndY = _hists2d[i]->GetBinContent(xBin, yBin);

        // Event detection prob for y
        double pY = _hists1d[_2dHistVars[i].second]->GetBinContent(
            _hists1d[_2dHistVars[i].second]->FindBin(point[_2dHistVars[i].second]));

        prob *= pXAndY / pY;
    }

    // Just make a cursory check that things are sensible
    if (prob > 1 || prob <= 0.0) {
        throw ChowLiu::BadEfficiency(prob, point);
    }

    return prob;
}
