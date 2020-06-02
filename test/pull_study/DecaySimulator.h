#ifndef DECAYSIMULATOR_HPP
#define DECAYSIMULATOR_HPP

#include <functional>
#include <random>
#include <utility>

#include "util.h"

/*
 * Class for running a Monte-Carlo simulaton of a D -> K3Pi experiment
 *
 * Calculate decay times for CF and DCS events with findCfDecayTimes() and findDcsDecayTimes()
 *
 * CF events generated according to to exponential decay
 * DCS events generated according to the (a + bt + ct^2) * exponential.
 *
 * Events selected from an exponential distribution
 *
 * As the DCS rate may not be bounded at large times, must provide a maximum time to the constructor.
 *
 * Based on the notes at this URL:
 *     http://www.columbia.edu/~ks20/4703-Sigman/4703-07-Notes-ARM.pdf
 */
class SimulatedDecays
{
  public:
    /*
     * Constructor
     * Seeds the random number generator used for generating times
     *
     * Sets the functions to be used as cf and dcs Rates; these can be e.g. lambdas
     *
     * Should also pass the function to be used to generate times (generateTime) and the PDF that this corresponds to
     * (generatingPDF). Note that NO CONSISTENCY CHECKS are performed between the time generation and the expected PDF;
     * it is entirely possible to pass in the "wrong" PDF for a given generator function. To test for this, a method
     * testGeneratingPDF is provided.
     *
     * To avoid seeding two separate random number generators which may not be indepdendent, a std::mt19937 should be
     * passed to this class which will be used to generate numbers for the accept-reject. This should be the same
     * mt19937 object used to create random numbers by generateTime().
     *
     * Also calculates the maximum ratio between our DCS decay and our model function; we need this to perform
     * accept-reject. Need to pass a std::pair of times to consider in order to calculate this. These should cover
     * the same range of times as generatingPDF.
     */
    SimulatedDecays(const std::function<double(void)> &  generateTime,
                    const std::function<double(double)> &generatingPDF,
                    const std::function<double(double)> &cfRate,
                    const std::function<double(double)> &dcsRate,
                    const std::pair<double, double> &    timeDomain,
                    const std::shared_ptr<std::mt19937> &generator);

    /*
     * Plot a graph of the generating PDF against a histogram of values generated using the provided time generation
     * function.
     *
     * Saves a plot to testPDF.pdf
     */
    void test(const size_t numPoints, const std::vector<double> &binLimits);

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
     * Max ratio between DCS rate and generating function
     * For UT
     */
    double maxDCSRatio(void);

    /*
     * Max ratio between CF rate and generating function
     * For UT
     */
    double maxCFRatio(void);

    /*
     * RS decay times
     */
    std::vector<double> RSDecayTimes{};

    /*
     * WS decay times
     */
    std::vector<double> WSDecayTimes{};

  private:
    /*
     * Generate a random number from 0 to 1 from a uniform distribution
     */
    double _getRandomUniform(void);

    /*
     * Find the maximum ratio of our rates to our exponential distribution and set _maxDCSRatio and _maxCFRatio
     */
    void _setMaxRatios(void);

    /*
     * Generate a random time according to the distribution passed in to the constructor
     */
    std::function<double(void)> _getRandomTime = [](void) {
        std::cerr << "generating fcn not set" << std::endl;
        throw D2K3PiException();
        return 0;
    };

    /*
     * The PDF we generate times from when _getRandomTime is called
     */
    std::function<double(double)> _generatingPDF = [](double x) {
        std::cerr << "generating pdf not set" << std::endl;
        throw D2K3PiException();
        return x;
    };

    /*
     * CF decay rate at a given time
     */
    std::function<double(double)> _cfRate = [](double x) {
        std::cerr << "cf rate not set" << std::endl;
        throw D2K3PiException();
        return x;
    };

    /*
     * DCS decay rate at a given time
     */
    std::function<double(double)> _dcsRate = [](double x) {
        std::cerr << "dcs rate not set" << std::endl;
        throw D2K3PiException();
        return x;
    };

    /*
     * Domain to consider
     */
    const double _minTime{0.0};
    const double _maxTime{0.0};

    /*
     * Maximum ratio of our DCS rate / exponential function
     * Needed to perform the acc-rej
     */
    double _maxDCSRatio{0.0};

    /*
     * Maximum ratio of our CF rate / exponential function
     * Needed to perform the acc-rej
     */
    double _maxCFRatio{0.0};

    /*
     * Mersenne Twister engine
     *
     * Initialised on construction
     */
    std::shared_ptr<std::mt19937> _gen = nullptr;

    /*
     * Uniform distribution to draw numbers from when _getRandomUniform() is called
     */
    std::uniform_real_distribution<double> _uniform;
};

#endif // DECAYSIMULATOR_HPP
