#include <iostream>

#include <TH1D.h>
#include <TH2D.h>

#include "ChowLiu.h"

namespace ChowLiu
{

HistogramProjections::HistogramProjections(const PhspBins& bins, const std::string& name)
    : _dimensionality(bins.size()), _bins(bins), _name(name)
{
    _1dhistograms = std::vector<std::unique_ptr<TH1D>>(_dimensionality);
    _2dhistograms = std::vector<std::unique_ptr<TH2D>>(_dimensionality * (_dimensionality - 1) / 2);
    _numBins      = std::vector<size_t>(_dimensionality);

    // Need to set the numbers of bins before creating histograms
    for (size_t i = 0; i < _dimensionality; ++i) {
        _numBins[i] = _bins[i].size() - 1;
    }

    for (size_t i = 0; i < _bins.size(); ++i) {
        std::string title1d = name + "1dhist" + std::to_string(i);
        _1dhistograms[i]    = std::make_unique<TH1D>(title1d.c_str(), title1d.c_str(), _numBins[i], _bins[i].data());

        for (size_t j = _bins.size() - 1; j > i; --j) {
            std::string title2d = name + "2dhist(" + std::to_string(i) + "," + std::to_string(j) + ")";
            _2dhistograms[_indexConversion(i, j)] = std::make_unique<TH2D>(
                title2d.c_str(), title2d.c_str(), _numBins[i], _bins[i].data(), _numBins[j], _bins[j].data());
        }
    }
}

void HistogramProjections::binPoint(const PhspPoint& point)
{
    if (point.size() != _dimensionality) {
        throw InvalidDimension();
    }

    for (size_t i = 0; i < _bins.size(); ++i) {
        _1dhistograms[i]->Fill(point[i]);

        for (size_t j = _bins.size() - 1; j > i; --j) {
            _2dhistograms[_indexConversion(i, j)]->Fill(point[i], point[j]);
        }
    }

    _numPoints++;
}

const HistogramProjections operator/(const HistogramProjections& numerator, const HistogramProjections& denominator)
{
    const PhspBins bins = numerator._bins;
    if (bins != denominator._bins) {
        throw BinMismatch();
    }
    std::string name = "ratio_" + numerator.getName() + "_" + denominator.getName();

    HistogramProjections result(bins, name);
    size_t               dimensionality = numerator.getDimensionality();

    for (size_t i = 0; i < dimensionality; ++i) {
        result._1dhistograms[i] = std::make_unique<TH1D>(numerator.get1dhistogram(i));
        TH1D denominator1dHist  = denominator.get1dhistogram(i);
        bool success            = result._1dhistograms[i]->Divide(&denominator1dHist);
        if (!success) {
            throw DivisionFailed();
        }

        for (size_t j = i + 1; j < dimensionality; ++j) {
            size_t index                = result._indexConversion(i, j);
            result._2dhistograms[index] = std::make_unique<TH2D>(numerator.get2dhistogram(i, j));
            TH2D denominator2dHist      = denominator.get2dhistogram(i, j);
            success                     = result._2dhistograms[index]->Divide(&denominator2dHist);
            if (!success) {
                throw DivisionFailed();
            }
        }
    }

    return result;
}

const TH1D HistogramProjections::get1dhistogram(const size_t i) const
{
    if (i > _dimensionality - 1) {
        throw HistogramNotFound();
    }
    return *_1dhistograms[i];
}

const TH2D HistogramProjections::get2dhistogram(const size_t i, const size_t j) const
{
    // Convert our pair of indices to the right 1d array index that we're using for storing our 2d histograms
    return *_2dhistograms[_indexConversion(i, j)];
}

size_t HistogramProjections::_indexConversion(const size_t i, const size_t j) const
{
    if (i == j) {
        std::cerr << "No 2d histogram of a variable against itself exists" << std::endl;
        throw HistogramNotFound();
    }

    // convert i, j -> smaller, larger
    size_t smaller = i < j ? i : j;
    size_t larger  = i > j ? i : j;

    if (larger > _dimensionality - 1) {
        throw HistogramNotFound();
    }

    // We store our 2d histograms in a 1d array
    // The array looks like {i,j} = {01, 02, 03 ..., 12, 13, 14... 23, 24 ...... (d-2)(d-1)}
    return (_dimensionality * (_dimensionality - 1) / 2) -
           (_dimensionality - smaller) * (_dimensionality - smaller - 1) / 2 + larger - smaller - 1;
}

Approximation::Approximation(const PhspBins& bins, const std::string& name)
    : HistogramProjections(bins, name), _graph(Graph(getDimensionality()))
{
    ;
}

void Approximation::makeApproximation(void)
{
    // Find the mutual information for each edge, add it to the graph
    for (size_t outNode = 0; outNode < getDimensionality(); ++outNode) {
        for (size_t inNode = outNode + 1; inNode < getDimensionality(); ++inNode) {
            TH2D   hist   = get2dhistogram(inNode, outNode);
            double weight = ChowLiu::mutual_info(&hist);
            _graph.addEdge(outNode, inNode, weight);
        }
    }

    // Find the directed MST and store it
    _directedTree = inTree(_root, _graph.getMaxSpanningTree());

    // Find what variables our histograms are in
    for (auto edges : _directedTree) {
        for (Edge edge : edges) {
            _2dHistVars.push_back(std::make_pair(edge.from(), edge.to()));
        }
    }

    // Populate our histograms
    for (size_t i = 0; i < getDimensionality(); ++i) {
        _hists1d.push_back(std::make_unique<TH1D>(get1dhistogram(i)));
    }

    for (auto pair : _2dHistVars) {
        // The HistogramProjections class only stores 2d histograms (i, j) for i<j
        // We may need to swap the axes of our histogram if we want (j, i) for our parametrisation
        TH2D hist = pair.first > pair.second ? swapAxes(get2dhistogram(pair.first, pair.second))
                                             : get2dhistogram(pair.first, pair.second);
        _hists2d.push_back(std::make_unique<TH2D>(hist));
    }

    _approximationMade = true;
}

double Approximation::value(const PhspPoint& point) const
{
    if (!_approximationMade) {
        throw ApproximationNotYetMade();
    }

    // Set the value to the 1d prob
    double prob1d = _hists1d[_root]->GetBinContent(_hists1d[_root]->FindBin(point[_root])) / (double)getNumPoints();
    double prob{prob1d};

    // Iterate over 2d histograms (there are d-1 of them), finding the conditional probabilities
    for (size_t i = 0; i < getDimensionality() - 1; ++i) {
        // The 2d histograms we stored are (x, y) for p(x|y)
        // The x and y bins our point is in on our 2d histogram
        size_t xBin = _hists2d[i]->GetXaxis()->FindBin(point[_2dHistVars[i].first]);
        size_t yBin = _hists2d[i]->GetYaxis()->FindBin(point[_2dHistVars[i].second]);

        // Number of events in both X and Y bin
        double nXY = _hists2d[i]->GetBinContent(xBin, yBin);

        // Number of events in the y bin
        // This could probably be streamlined a bit with yBin
        double nY = _hists1d[_2dHistVars[i].second]->GetBinContent(
            _hists1d[_2dHistVars[i].second]->FindBin(point[_2dHistVars[i].second]));

        prob *= nXY / nY;
    }

    // Just make a cursory check that things are sensible
    if (prob > 1 || prob <= 0.0) {
        throw BadEfficiency(prob, point);
    }

    return prob;
}

double entropy(const TH1D* const hist)
{
    size_t totalEntries = hist->Integral();
    double ent{0.0}; // entropy

    for (size_t bin = 1; bin <= (size_t)hist->GetNbinsX(); ++bin) { // ROOT bin indexing starts at 1
        double binContent = hist->GetBinContent(bin);
        if (binContent > FLT_EPSILON) {
            double prob = binContent / totalEntries;
            ent -= prob * std::log(prob);
        }
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

    size_t numEvents = xHist->Integral();
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

TH2D swapAxes(const TH2D& other)
{
    size_t nBinsX = other.GetNbinsX();
    size_t nBinsY = other.GetNbinsY();

    const double*     xBins = other.GetXaxis()->GetXbins()->GetArray();
    const double*     yBins = other.GetYaxis()->GetXbins()->GetArray();
    std::string       oldName(other.GetName());
    const std::string name("swap_" + oldName);

    TH2D swappedHist = TH2D(name.c_str(), name.c_str(), nBinsY, yBins, nBinsX, xBins);

    for (size_t i = 1; i <= nBinsX; ++i) {
        for (size_t j = 1; j <= nBinsY; ++j) {
            swappedHist.SetBinContent(i, j, other.GetBinContent(j, i));
        }
    }
    return swappedHist;
}

} // namespace ChowLiu
