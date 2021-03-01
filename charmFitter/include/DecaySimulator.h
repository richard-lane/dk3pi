#ifndef DECAYSIMULATOR_HPP
#define DECAYSIMULATOR_HPP

#include <random>

#include "fitterUtil.h"

/*
 * Class for running a Monte-Carlo simulaton of a D -> K3Pi experiment
 *
 * Calculate decay times for RS and WS events with findCfDecayTimes() and findDcsDecayTimes()
 *
 * RS events generated according to to exponential decay
 * WS events generated according to the (a + bt + ct^2) * exponential.
 *
 * WS events selected using accept-reject, generated from an exponential distribution.
 * Based on the notes at this URL:
 *     http://www.columbia.edu/~ks20/4703-Sigman/4703-07-Notes-ARM.pdf
 */
class SimulatedDecays
{
  public:
    /*
     * As the DCS rate may not be bounded at large times (at least, when we use our parametrisation that is only valid
     * at small times...), the user must provide a maximum time to the constructor.
     *
     * Must also provide an initialised random number generator
     */
    SimulatedDecays(const std::pair<double, double>& timeDomain,
                    const FitterUtil::DecayParams_t& decayParams,
                    std::mt19937& rng); // I've read in places that this isn't a good solution... but what is?

    /*
     * Generate a decay time from the RS distribution
     */
    double rsPoint();

    /*
     * Generate a decay time from the WS distribution
     */
    double wsPoint();

    /*
     * Generate a collection of decay times from the RS distribution
     */
    std::vector<double> rsDecayTimes(const size_t numPoints);

    /*
     * Generate a collection of decay times from the WS distribution
     */
    std::vector<double> wsDecayTimes(const size_t numPoints);

  private:
    /*
     * The decay rates at a given time
     *
     * These are not scaled: the RS rate is ~ (e^-kt) and the WS ~ (a + bt + ct^2)(e^-kt).
     * To produce a sensible ratio of RS/WS points, the user should calculate how many RS and WS events they want
     */
    double _rsRate(const double time) const;
    double _wsRate(const double time) const;

    /*
     * Find the maximum ratio of our rates
     *
     * Differentiating a + bt + ct^2, we find a turning point at t = -b/2c.
     * This means our maximum ratio is either at _minTime, _maxTime or -b/2c if it is in range.
     * We just need to check the ratio at all of these times to find where it is maximal.
     */
    double _findMaxRatio(void) const;

    /*
     * Check whether a time is accepted as being from the WS distribution
     */
    bool _isAccepted(const double time);

    /*
     * Time domain to consider
     */
    const double _minTime{0.0};
    const double _maxTime{0.0};

    /*
     * Parameters describing these decays
     */
    FitterUtil::DecayParams_t _decayParams{};

    /*
     * RNG
     * Should probably be a reference instead of a pointer
     */
    std::mt19937* _gen{nullptr};

    /*
     * Max ratios of our WS pdfs to our exponential distribution
     * Needed to perform the accept-reject
     */
    double _maxWSRatio{0.0};

    /*
     * Uniform distribution to draw numbers from when performing the accept-reject
     *
     * Also used to generate exponentially-distributed times
     */
    std::uniform_real_distribution<double> _uniform{};
};

#endif // DECAYSIMULATOR_HPP
