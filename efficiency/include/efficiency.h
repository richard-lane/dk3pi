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
typedef std::array<std::vector<double>, 5> PhspBins;

/*
 * Parametrise a phase space point by its 5d coordinates (p0, p1, p2, p3, p4)
 */
typedef std::array<double, 5> PhspPoint;

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
    const TH1D get1dhistogram(const unsigned short i) const;

    /*
     * Get the 2d histogram relating the i and j'th parameters
     *
     * i must not equal j
     */
    const TH2D get2dhistogram(const unsigned short i, const unsigned short j) const;

    /*
     * Bin a point
     */
    void binPoint(const PhspPoint& point);

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
    // Hard-coded array-lengths here are probably fine since my analysis only cares about 5d efficiencies...
    /*
     * Histograms representing {p0, p1, p2, p3, p4}
     */
    std::array<std::unique_ptr<TH1D>, 5> _1dhistograms{};

    /*
     * Histograms representing {p0vp1, p0vp2, p0vp3, p0vp4, p1vp2, p1vp3 ...}
     */
    std::array<std::unique_ptr<TH2D>, 10> _2dhistograms{};

    /*
     * Number of bins along each of our {p0, p1, p2, p3, p4} axes
     */
    std::array<size_t, 5> _numBins{};

    /*
     * Array of 5 vectors of bins to use
     */
    PhspBins _bins{};

    /*
     * Convert an (i, j) index pair to an index between 0 and 9 for our 2d histogram array
     */
    unsigned short _indexConversion(const unsigned short i, const unsigned short j) const;
};

/*
 * Class for estimating the phase space dependent detection efficiency of a process given truth- and detector-level
 * events
 */
class Efficiency
{
  public:
    explicit Efficiency(const std::vector<dDecay_t>&            detectedEvents,
                        const std::vector<dDecay_t>&            generatedEvents,
                        const std::vector<std::vector<double>>& bins);

    /*
     * Return the value of the efficiency function at a phase space point (m12, m23, m34, m123, m234)
     */
    double value(const std::vector<double>& invMasses) const;

  private:
    std::vector<dDecay_t>            _detectedEvents;
    std::vector<dDecay_t>            _generatedEvents;
    std::vector<std::vector<double>> _bins;

    // It might be better to have 10 vectors of m12 values etc.
    std::vector<std::vector<double>> _detectedInvMasses;
    std::vector<std::vector<double>> _generatedInvMasses;

    /*
     * Find the invariant masses of each event
     */
    void _findInvariantMasses(void);
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
