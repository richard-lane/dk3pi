#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include <memory>
#include <vector>

#include <TH1D.h>
#include <TH2D.h>

#include "ChowLiu.h"
#include "efficiencyUtil.h"
#include "graphTheory.h"
#include "util.h"

/*
 * Class for estimating the phase space dependent detection efficiency of a process given truth- and detector-level
 * events
 *
 * Give it a load of truth- and detector-level events, call efficiencyParametrisation() once all the data is added and
 * then call value() to find the value of the efficiency at a phase space point (or rather, the value of the efficiency
 * in the bin where that phsp point lives)
 */
class ChowLiuEfficiency
{
  public:
    /*
     * Tell us the bins we're using and the root node to use for the approximation
     */
    explicit ChowLiuEfficiency(const PhspBins& bins, const size_t root = 0);

    /*
     * Return the value of the efficiency function at a phase space point
     */
    double value(const PhspPoint& point) const;

    /*
     * Add a monte carlo event
     */
    void addMCEvent(const PhspPoint& point);

    /*
     * Add a generated "truth-level" event
     */
    void addGeneratedEvent(const PhspPoint& point);

    /*
     * Add a collection of MC events
     */
    template <typename ContainerType> void addMCEvents(const ContainerType& points)
    {
        for (auto point : points) {
            addMCEvent(point);
        }
    }

    /*
     * Add a collection of truth level events
     */
    template <typename ContainerType> void addGeneratedEvents(const ContainerType& points)
    {
        for (auto point : points) {
            addGeneratedEvent(point);
        }
    }

    /*
     * Work out the efficiency projections in 1/2d from our MC and truth-level datasets, find the best approximation to
     * the efficiency function in 5d, and set up all the bits so that this->value() works
     */
    void efficiencyParametrisation(void);

  private:
    PhspBins                                       _bins;
    std::unique_ptr<ChowLiu::HistogramProjections> _detectedEvents  = nullptr;
    std::unique_ptr<ChowLiu::HistogramProjections> _generatedEvents = nullptr;

    /*
     * The ratio _detectedEvents / _generatedEvents; i.e. the efficiency projections
     */
    ChowLiu::HistogramProjections _ratio = ChowLiu::HistogramProjections(_bins, "none");

    /*
     * The node to use as root of our approximation
     */
    size_t _root{0};

    /*
     * Flag tracking whether our efficiency parametrisation has been made yet
     */
    bool _approximationMade{false};

    /*
     * Dimensionality of the probability distribution we are approximating
     */
    size_t _dimensionality{0};

    /*
     * Graph object used to calculate which variables we want to keep
     */
    std::unique_ptr<Graph> _graph = nullptr;

    /*
     * Adjacency list representation of the (directed) graph used to calculate which pairs of variables to keep
     */
    std::vector<std::list<Edge>> _directedTree;

    /*
     * The 1d histograms we use in our efficiency approximation
     */
    std::vector<std::unique_ptr<TH1D>> _hists1d{};

    /*
     * The 2d histograms we'll use in our efficiency approximation
     * Should be a histogram of i vs j if we need a term like p(i|j)
     */
    std::vector<std::unique_ptr<TH2D>> _hists2d{};

    /*
     * Which variables our 2d histograms are in
     * e.g. if we are parametrising as e(0|1)e(2|1)e(3|2), this should be {(0, 1), (2, 1), (3, 2)}
     */
    std::vector<std::pair<size_t, size_t>> _2dHistVars{0};
};

#endif // EFFICIENCY_H
