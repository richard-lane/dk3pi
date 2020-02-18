/*
 * Class for simulating the D->K3Pi decays according to time-dependent 2nd order polynomials
 */
#ifndef DECAYSIMULATOR_HPP
#define DECAYSIMULATOR_HPP

#include <utility>

#include "MCGenerator.h"

/*
 * Struct encapsulating the parameters needed to simulate decays.
 * Mixing params x, y
 * Phase space params r, Im(Z) and Re(Z)
 * Particle data \Gamma for the D meson
 *
 * Parameters are all initalised to zero by defauly and should be set by the user.
 */
typedef struct DecayParameters {
    // Mixing params
    double x{0.0};
    double y{0.0};

    // Phase space params
    double r{0.0};
    double z_im{0.0};
    double z_re{0.0};

    // Particle data
    double width{0.0};
} DecayParams_t;

class SimulatedDecays : public MCGenerator
{
  public:
    /*
     * Set the allowed time and decay rate values for our simulated decay.
     */
    SimulatedDecays(const std::pair<double, double> &timeRange,
                    const std::pair<double, double> &decayRateRange,
                    const DecayParams_t &            DecayParams);

    /*
     * Check whether a point is accepted as being from one of the distributions.
     * Boolean flag rightSign indicates whether to use the RS or WS decay model.
     */
    bool isAccepted(const double xVal, const double yVal, bool rightSign);

    /*
     * Generate a vector of numEvents decay times representing DCS decay times.
     *
     * Sets WSDecayTimes
     */
    void findDcsDecayTimes(size_t numEvents);

    /*
     * Generate a vector of numEvents decay times representing CF decay times.
     *
     * Sets RSDecayTimes
     */
    void findCfDecayTimes(size_t numEvents);

    /*
     * Plot graphs of the WS and RS decay rates
     */
    void plotRates(void);

    /*
     * RS decay times
     */
    std::vector<double> RSDecayTimes{};

    /*
     * WS decay times
     */
    std::vector<double> WSDecayTimes{};

    /*
     * Find the RS decay rate at a given time.
     * We only care about the ratio of RS to WS decays, so an overall factor (of B^2 as defined in papers) is omitted.
     * Public for UT
     */
    double _rightSignDecayRate(const double time);

    /*
     * Find the WS decay rate at a given time.
     * We only care about the ratio of RS to WS decays, so an overall factor (of B^2 as defined in papers) is omitted.
     * Public for UT
     */
    double _wrongSignDecayRate(const double time);

  private:
    /*
     * Mixing parameters
     */
    DecayParams_t _DecayParams;
};

#endif // DECAYSIMULATOR_HPP
