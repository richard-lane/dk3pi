#ifndef CHOWLIU_H
#define CHOWLIU_H

#include <memory>
#include <string>
#include <vector>

#include <TH2D.h>

#include "efficiencyUtil.h"
#include "graphTheory.h"

namespace ChowLiu
{

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
 * Class for binning points in multiple 1- and 2-d histograms
 *
 */
class HistogramProjections
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
    explicit HistogramProjections(const PhspBins& bins, const std::string& name);

    /*
     * Bin a point
     */
    void binPoint(const PhspPoint& point);

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

    size_t getNumPoints(void) const { return _numPoints; };

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

  protected:
    /*
     * Number of dimensions of our phsp
     */
    size_t _dimensionality{0};

    /*
     * Number of bins along each of our {p0, p1, p2, ...} axes
     */
    std::vector<size_t> _numBins{};

    /*
     * Bins that we're using
     */
    PhspBins _bins{};

  private:
    /*
     * 1d Histograms
     */
    std::vector<std::unique_ptr<TH1D>> _1dhistograms{};

    /*
     * Histograms representing {p0vp1, p0vp2, p0vp3, p0vp4, p1vp2, p1vp3 ...}
     *
     * Indexed according to _indexConversion(i, j); (pi vs pj) found at _indexConversion(i, j)
     */
    std::vector<std::unique_ptr<TH2D>> _2dhistograms{};

    /*
     * Something to label this object with
     */
    std::string _name;

    /*
     * I just don't trust ROOT's histograms to correctly track how many points I've entered
     */
    size_t _numPoints{0};

    /*
     * Convert an (i, j) index pair to an index for our array of 2d histograms
     */
    size_t _indexConversion(const size_t i, const size_t j) const;
};

/*
 * Ratio of two projections
 */
const HistogramProjections operator/(const HistogramProjections& numerator, const HistogramProjections& denominator);

struct BinMismatch : public std::exception {
    const char* what() const throw() { return "Numerator and denominator bins do not match"; }
};

struct ApproximationNotYetMade : public std::exception {
    const char* what() const throw()
    {
        return "Cannot retrieve efficiency value; first call this->efficiencyParametrisation";
    }
};

struct BadProbability : public std::exception {
    BadProbability(const double value, const PhspPoint& point)
    {
        _msg = "Probability " + std::to_string(value) + " returned at point (";
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
 * Class for making the Chow-Liu approximation for a probability distribution
 *
 * Our distribution is represented by histograms- the distribution we're interested in here is p(x) for an event to
 * occur at x
 *
 * Usage:
 * Add points to the histograms using this->binPoint() or this->binPoints(), then make the approximation using
 * this->makeApproximation An approximation to the probability at p can then be found from this->value(p)
 */
class Approximation : public HistogramProjections
{
  public:
    /*
     * The bins to use + the name to label this dataset by
     */
    Approximation(const PhspBins& bins, const std::string& name);

    /*
     * Work out which variables to use in our Chow-Liu approximation and set up the bits so that this->value() works
     */
    void makeApproximation(void);

    /*
     * Return the probability at a point
     */
    double value(const PhspPoint& point) const;

  private:
    /*
     * The node to use as root of our approximation
     */
    size_t _root{0};

    /*
     * Flag tracking whether our efficiency parametrisation has been made yet
     */
    bool _approximationMade{false};

    /*
     * Graph object used to calculate which variables we want to keep
     */
    Graph _graph;

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

} // namespace ChowLiu

#endif // CHOWLIU_H
