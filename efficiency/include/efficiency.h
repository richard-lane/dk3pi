#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include <array>
#include <memory>
#include <vector>

#include <TH1D.h>
#include <TH2D.h>

#include "efficiencyUtil.h"
#include "util.h"

/*
 * Bins that we will use for our efficiency parametrisation
 */
typedef std::vector<std::vector<double>> PhspBins;

/*
 * Parametrise a phase space point by its coordinates (p0, p1, p2, p3, p4)
 */
typedef std::vector<double> PhspPoint;

struct InvalidDimension : public std::exception {
    const char* what() const throw() { return "Dimension of point does not match dimension of phase space"; }
};

struct HistogramNotFound : public std::exception {
    const char* what() const throw() { return "No Histogram exists at the specified index"; }
};

/*
 * Class for binning phase space points in 5 1d histograms
 * Bin edges and phase space points should be of type T
 *
 * They'll probably all be doubles but templates are just more fun, aren't they
 *
 * Must provide PhspBins to the ctor
 *
 */
class EfficiencyBinning
{
  public:
    /*
     * Initialise our bins but don't put anything in them yet
     *
     * Each vector of bin limits should contain the low edge of each bin, plus the high edge of the last bin
     */
    explicit EfficiencyBinning(const PhspBins& bins);

    /*
     * Initialise our bins and put a single point in them
     */
    EfficiencyBinning(const PhspBins& bins, const PhspPoint& point);

    /*
     * Bin a point
     */
    void binPoint(const PhspPoint& point);

    /*
     * Initialise our bins and put several points in them
     */
    template <typename ContainerType>
    EfficiencyBinning(const PhspBins& bins, const ContainerType& points) : EfficiencyBinning(bins)
    {
        binPoints(points);
    }

    /*
     * Get the i'th 1d histogram
     *
     * Indexed as {0:p0, 1:p1, ... }
     */
    const TH1D get1dhistogram(const size_t i) const;

    /*
     * Get the 2d histogram relating the i and j'th parameters
     *
     * i must not equal j
     */
    const TH2D get2dhistogram(const size_t i, const size_t j) const;

    /*
     * Bin a collection of points
     */
    template <typename ContainerType> void binPoints(const ContainerType& points)
    {
        {
            for (const auto& point : points) {
                binPoint(point);
            }
        }
    }

  private:
    /*
     * Number of dimensions of our phsp
     */
    size_t _dimensionality{0};

    /*
     * Histograms representing {p0, p1, p2, p3, p4}
     */
    std::vector<std::unique_ptr<TH1D>> _1dhistograms{};

    /*
     * Histograms representing {p0vp1, p0vp2, p0vp3, p0vp4, p1vp2, p1vp3 ...}
     * Indexed in a weird order; histogram (pi vs pj) found at _2dhistograms[_indexConversion(i, j)]
     */
    std::vector<std::unique_ptr<TH2D>> _2dhistograms{};

    /*
     * Number of bins along each of our {p0, p1, p2, ...} axes
     */
    std::vector<size_t> _numBins{};

    /*
     * Bins that we're using
     */
    PhspBins _bins{};

    /*
     * Convert an (i, j) index pair to an index for our 2d histogram array
     */
    size_t _indexConversion(const size_t i, const size_t j) const;
};

/*
 * Class for estimating the phase space dependent detection efficiency of a process given truth- and detector-level
 * events
 */
class ChowLiuEfficiency
{
  public:
    explicit ChowLiuEfficiency(const PhspBins& bins);

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
    EfficiencyBinning _detectedEvents;
    EfficiencyBinning _generatedEvents;
    PhspBins          _bins;

    /*
     * Flag tracking whether our efficiency parametrisation has been made yet
     */
    bool _approximationMade{false};

    /*
     * The 1d histogram we use in our efficiency approximation
     */
    std::unique_ptr<TH1D> _hist1d = nullptr;

    /*
     * Which variable our 1d histogram is in
     */
    size_t _1dHistVar{0};

    /*
     * The 2d histograms we'll use in our efficiency approximation
     */
    std::vector<std::unique_ptr<TH2D>> _hists2d{};

    /*
     * Which variables our 2d histograms are in
     * e.g. if we are parametrising as e(0|1)e(2|1)e(3|2), this should be {(0, 1), (2, 1), (3, 2)}
     */
    std::vector<std::pair<size_t, size_t>> _2dHistVars{0};
};

/*
 * Work out the entropy associated with a histogram
 *
 * entropy = Sum(-P logP) for probability of bin contents P
 *
 * Uses log base e, so entropy is returned in natural units
 */
double entropy(const TH1D* const hist);

/*
 * Calculate the mutual information shared by two labellings of the same data
 *
 * i.e. bin a histogram in 2d according to the variables (X1, X2) and pass it to this function to find the mutual
 * information between X1 and X2
 *
 * Takes a pointer because thats what ROOT likes to do
 */
double mutual_info(const TH2D* const histogram2d);

#endif // EFFICIENCY_H
