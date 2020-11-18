#include "scriptUtils.h"

#include <TH1D.h>

HistogramSlices::HistogramSlices(const std::string&               title,
                                 const size_t                     numSlices,
                                 const size_t                     numBins,
                                 const std::pair<double, double>& histLimits,
                                 const std::pair<double, double>& sliceLimits,
                                 const size_t                     plotVarIndex,
                                 const size_t                     sliceVarIndex)
    : _sliceHist(TH1D((title + "Slice Hist").c_str(), "Slice Hist", numSlices, sliceLimits.first, sliceLimits.second)),
      _plotVarIndex(plotVarIndex), _sliceVarIndex(sliceVarIndex)
{
    _slices.reserve(numSlices);
    for (size_t i = 0; i < numSlices; ++i) {
        const std::string thisTitle = "Slice " + std::to_string(i) + ": " + title;
        _slices.push_back(TH1D(thisTitle.c_str(), thisTitle.c_str(), numBins, histLimits.first, histLimits.second));
        _slices[i].SetStats(false);
    }
}

void HistogramSlices::add(const PhspPoint& point, const double wt)
{
    // Find the right histogram to put it in
    // ROOT bin indexing starts from 1
    const size_t slice = _sliceHist.FindBin(point[_sliceVarIndex]) - 1;

    // Add this point to it
    _slices[slice].Fill(point[_plotVarIndex], wt);

    _numPoints += wt;
}

void HistogramSlices::add(const std::vector<PhspPoint>& points, const std::vector<double>* wts)
{
    if (wts) {
        assert(points.size() == wts->size());
        for (size_t i = 0; i < points.size(); ++i) {
            this->add(points[i], (*wts)[i]);
        }

    } else {
        for (const auto& point : points) {
            this->add(point);
        }
    }
}

void HistogramSlices::setColour(const EColor colour)
{
    for (auto& hist : _slices) {
        hist.SetLineColor(colour);
    }
}

void plotSlices(const std::string&               path,
                std::vector<HistogramSlices>&    slices,
                const std::vector<std::string>&  plotOptions,
                const std::vector<std::string>&  labels,
                const std::pair<double, double>& range)
{
    assert(slices.size() == plotOptions.size());
    assert(slices.size() == labels.size());

    size_t numSlices = slices[0]._slices.size();
    for (const auto& histSlice : slices) {
        assert(numSlices == histSlice._slices.size());
    }

    for (size_t i = 0; i < numSlices; ++i) {
        const std::string    slicePath{path + std::to_string(i) + ".png"};
        std::vector<TH1D*>   hists{};
        util::LegendParams_t legend{0.7, 0.9, 0.7, 0.9};

        for (auto& histSlice : slices) {
            // Scale hists so they all look nice
            histSlice._slices[i].Scale(1 / histSlice._numPoints);
            hists.push_back(&histSlice._slices[i]);

            // Set Y axis range
            histSlice._slices[i].GetYaxis()->SetRangeUser(range.first, range.second);
        }
        // Plot hists
        util::saveObjectsToFile<TH1D>(
            std::vector<TObject*>(hists.begin(), hists.end()), plotOptions, labels, slicePath, legend);
    }

    // Rescale. Just in case we want to plot again after this fcn returns
    // If there is a failure before now this will leave the hist slices in a broken state, but that isn't likely to be
    // important
    for (size_t i = 0; i < numSlices; ++i) {
        for (auto& histSlice : slices) {
            histSlice._slices[i].Scale(1 / histSlice._numPoints);
        }
    }
}

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