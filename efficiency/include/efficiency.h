#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include <vector>

#include <TH2D.h>

#include "efficiencyUtil.h"

/*
 * Bins that we will use for our efficiency parametrisation
 *
 * Made it a template cus why not
 */
template <typename T> struct PhspBins {
    std::vector<T> m12_bins;
    std::vector<T> m23_bins;
    std::vector<T> m34_bins;
    std::vector<T> m123_bins;
    std::vector<T> m234_bins;
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
