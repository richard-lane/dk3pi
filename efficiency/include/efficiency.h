#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include <array>
#include <memory>
#include <vector>

#include <TH1D.h>
#include <TH2D.h>

#include "efficiencyUtil.h"
#include "graphTheory.h"
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

struct DivisionFailed : public std::exception {
    const char* what() const throw() { return "For ROOT reasons, division failed"; }
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
     *
     * Pass in a name to label this object, because two histograms existing with the same name causes ROOT issues
     * possibly
     */
    explicit EfficiencyBinning(const PhspBins& bins, const std::string& name);

    /*
     * Initialise our bins and put a single point in them
     */
    EfficiencyBinning(const PhspBins& bins, const PhspPoint& point, const std::string& name);

    /*
     * Bin a point
     */
    void binPoint(const PhspPoint& point);

    /*
     * Initialise our bins and put several points in them
     */
    template <typename ContainerType>
    EfficiencyBinning(const PhspBins& bins, const ContainerType& points, const std::string& name)
        : EfficiencyBinning(bins, name)
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

    size_t getDimensionality(void) const { return _dimensionality; };

    const std::string getName(void) const { return _name; };

    /*
     * Take the ratio of two EfficiencyBinning instances
     *
     * Quite hacky and weird; uses the underlying histograms directly then sets them
     */
    friend const EfficiencyBinning operator/(const EfficiencyBinning& numerator, const EfficiencyBinning& denominator);

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
     * Something to label this object with
     */
    std::string _name;

    /*
     * Convert an (i, j) index pair to an index for our 2d histogram array
     */
    size_t _indexConversion(const size_t i, const size_t j) const;
};

/*
 * Ratio of two projections
 */
const EfficiencyBinning operator/(const EfficiencyBinning& numerator, const EfficiencyBinning& denominator);

struct BinMismatch : public std::exception {
    const char* what() const throw() { return "Numerator and denominator bins do not match"; }
};

struct ApproximationNotYetMade : public std::exception {
    const char* what() const throw()
    {
        return "Cannot retrieve efficiency value; first call this->efficiencyParametrisation";
    }
};

struct BadEfficiency : public std::exception {
    BadEfficiency(const double value, const PhspPoint& point)
    {
        _msg = "Efficiency value " + std::to_string(value) + " returned at point (";
        for (auto coord : point) {
            _msg += std::to_string(coord) + ", ";
        }
        // Remove the unnecessary ", " from the end of the last coord
        _msg.erase(_msg.length() - 2);
        _msg += ")";
    }
    const char* what() const throw() { return _msg.c_str(); }

  private:
    std::string _msg; // Probably bad since constructing a string could throw. but idc
};

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
    PhspBins                           _bins;
    std::unique_ptr<EfficiencyBinning> _detectedEvents  = nullptr;
    std::unique_ptr<EfficiencyBinning> _generatedEvents = nullptr;

    /*
     * The ratio _detectedEvents / _generatedEvents; i.e. the efficiency projections
     */
    EfficiencyBinning _ratio = EfficiencyBinning(_bins, "none");

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

/*
 * Take a 2d histogram and swap x<->y axes
 */
TH2D swapAxes(const TH2D& other);

#endif // EFFICIENCY_H
