#ifndef EFFICIENCY_SCRIPT_UTILS_H
#define EFFICIENCY_SCRIPT_UTILS_H

#include <iostream>
#include <memory>
#include <random>
#include <string>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLegend.h>

#include "ReadRoot.h"
#include "efficiency.h"
#include "efficiencyUtil.h"
#include "util.h"

struct EventDetectionProbNotNormalised : public std::exception {
    EventDetectionProbNotNormalised(const double prob)
        : _msg("Probability " + std::to_string(prob) + " not between 0 and 1")
    {
        ;
    }

    const char* what() const throw() { return _msg.c_str(); }

  private:
    const std::string _msg;
};

/*
 * Helper for representing a multidimensional distribution as slices through a histogram
 *
 */
class HistogramSlices
{
  public:
    /*
     * Hists named <title>
     *
     * Titles should follow ROOT's convention of the title string being semicolon-separated:
     *     "<global plot title>;<x axis title>;<y axis title>"
     */
    HistogramSlices(const std::string&               title,
                    const size_t                     numSlices,
                    const size_t                     numBins,
                    const std::pair<double, double>& histLimits,
                    const std::pair<double, double>& sliceLimits,
                    const size_t                     plotVarIndex,
                    const size_t                     sliceVarIndex);

    /*
     * Add a point to the histograms
     *
     * Weights default to 1 if unspecified
     */
    void add(const PhspPoint& point, const double wt = 1);
    void add(const std::vector<PhspPoint>& points, const std::vector<double>* weights = nullptr);

    void setColour(const EColor colour);

  private:
    /*
     * Plots slices on the same canvas to <path>0.png, <oath>1.png, <path>2.png,...
     *
     * Have to set Y limits here, as if you set limits then scale it breaks.
     *
     * yep
     */
    friend void plotSlices(const std::string&               path,
                           std::vector<HistogramSlices>&    slices,
                           const std::vector<std::string>&  plotOptions,
                           const std::vector<std::string>&  labels,
                           const std::pair<double, double>& range);

    std::vector<TH1D> _slices{};

    /*
     * Hist for tracking which slice a point belongs in
     */
    TH1D _sliceHist{};

    /*
     * Indices of variables we're taking the slices across and we're making the histograms in
     */
    size_t _plotVarIndex{};
    size_t _sliceVarIndex{};

    /*
     * Floating point since the points can have weights
     */
    double _numPoints{0};
};

/*
 * Provide a random number generator, a vector of D decay events and a function that provides an event's detection
 * probability
 *
 * The probability of detecting an event should be a number between 0 and 1
 *
 * Removes the undetected events from the vector in-place.
 */
void applyEfficiency(std::mt19937* const                    generator,
                     const std::function<double(dDecay_t)>& eventDetectionProb,
                     std::vector<dDecay_t>&                 events);

/*
 * Return the 5d parametrisation of a phase space point
 */
PhspPoint parametrisation(const dDecay_t& decay);

inline double simpleEfficiency([[maybe_unused]] const dDecay_t& event)
{
    return 0.5;
}

/*
 * A nice efficiency that gets recovered well
 */
inline double niceEfficiency(const dDecay_t& event)
{
    return invariantMass({event.kParams, event.pi1Params}) / 3;
}

/*
 * A less-nice efficiency that doesn't get recovered as nicely
 */
inline double awkwardEfficiency(const dDecay_t& event)
{
    std::vector<double> params = parametrisation(event);
    double              e{1};
    for (int i = 0; i < 3; ++i) {
        e *= params[i] / 2;
    }
    return 5 * e;
}

/*
 * An efficiency on the total pT of the k
 */
inline double pTEfficiency(const dDecay_t& event)
{
    return std::sqrt(pT(event.kParams));
}

/*
 * Return some sensible looking bin limits
 */
PhspBins findBins(void);

#endif // EFFICIENCY_SCRIPT_UTILS_H
