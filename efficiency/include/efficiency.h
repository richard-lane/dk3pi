#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include <vector>

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
    explicit Efficiency(const std::vector<dDecay_t>& detectedEvents,
                        const std::vector<dDecay_t>& generatedEvents,
                        const PhspBins<double>&      bins);

    /*
     * Return the value of the efficiency function at the phase space point described by an event
     */
    double value(const dDecay_t& event) const;

  private:
    std::vector<dDecay_t> _detectedEvents;
    std::vector<dDecay_t> _generatedEvents;
    PhspBins<double>      _bins;

    // It might be better to have 10 vectors of m12 values etc.
    std::vector<std::vector<double>> _detectedInvMasses;
    std::vector<std::vector<double>> _generatedInvMasses;

    /*
     * Find the invariant masses of each event
     */
    void _findInvariantMasses(void);
};

#endif // EFFICIENCY_H
