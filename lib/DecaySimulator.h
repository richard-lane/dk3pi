/*
 * Class for simulating the D->K3Pi decays according to time-dependent 2nd order polynomials
 */
#ifndef DECAYSIMULATOR_HPP
#define DECAYSIMULATOR_HPP

#include <random>
#include <utility>

#include "util.h"

/*
 * Class for running a Monte-Carlo simulaton of a D -> K3Pi experiment
 * Calculate decay times for CF and DCS events with findCfDecayTimes() and findDcsDecayTimes()
 *
 * CF events generated according to exponential decay; DCS events generated according to the (a + bt + ct^2) *
 * exponential.
 *
 * As the DCS rate may not be bounded at large times, must provide a maximum time to the constructor.
 */
class SimulatedDecays
{
  public:
    /*
     * Set the allowed time and decay rate values for our simulated decay, and seed our random number generator.
     *
     * Also calculate the maximum ratio between our DCS decay and our model function; we need this to perform
     * accept-reject
     */
    SimulatedDecays(const double maxTime, const DecayParams_t &DecayParams);

    /*
     * Check whether a time is accepted as being from one of the distributions, using a number selected froma uniform
     * distribution
     *
     * Boolean flag rightSign indicates whether to use the RS or WS decay model.
     */
    bool isAccepted(const double time, const double uniformVal, bool rightSign);

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
     * Plot histograms of the number of decay rates in each time bin
     */
    void plotRates(const std::vector<double> &binLimits);

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
     *
     * Public for UT
     *
     * TODO remove
     */
    double rightSignDecayRate(const double time);

    /*
     * Find the WS decay rate at a given time.
     * We only care about the ratio of RS to WS decays, so an overall factor (of B^2 as defined in papers) is omitted.
     *
     * Public for UT
     *
     * TODO remove
     */
    double wrongSignDecayRate(const double time);

  private:
    /*
     * Generate a random number from 0 to 1 from a uniform distribution
     */
    double _getRandomUniform(void);

    /*
     * Generate a random time, exponentially distributed with characteristic time 1/_DecayParams.width, between 0 and
     * _maxTime
     */
    double _getRandomTime(void);

    /*
     * Find the maximum ratio of our DCS ratio to our exponential distribution and set _maxDCSRatio
     *
     * Parameter c in (a + bt + ct^2) is necessarily positive (c = 0.25 * x^2 * y^2 * width*2), so the maximum will
     * either be at t=0 or t=maxTime
     */
    void _setMaxDCSRatio(void);

    /*
     * Times to generate up to
     * Cannot be arbitrarily large as our (a + bt + ct^2) DCS/CF rate approximation is only valid at low times, and
     * hence is not bounded at large times
     */
    double _maxTime{0.0};

    /*
     * Maximum ratio of our DCS rate / exponential function
     * Needed to perform the acc-rej
     */
    double _maxDCSRatio{0.0};

    /*
     * Mersenne Twister engine
     *
     * Initialised on construction
     */
    std::mt19937 _gen;

    /*
     * Uniform distribution to draw numbers from when _getRandomUniform() is called
     */
    std::uniform_real_distribution<double> _uniform;

    /*
     * Mixing parameters
     */
    DecayParams_t _DecayParams;
};

#endif // DECAYSIMULATOR_HPP
