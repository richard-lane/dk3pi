/*
 * Monte-Carlo generation of points in a distribution
 */

#ifndef MCGENERATOR_HPP
#define MCGENERATOR_HPP

#include <random>

/*
 * Base class for generation of points using Monte Carlo
 *
 * Contains methods for generating from a distribution and running accept/reject
 * Intented to be a parent class to useful implementations
 */
class MCGenerator
{
  public:
    MCGenerator(const double XRangeMin = 0,
                const double XRangeMax = 1,
                const double YRangeMin = 0,
                const double YRangeMax = 1);

    /*
     * Set the minimum X value that can be generated
     */
    void _setMinXValue(const double value);

    /*
     * Set the maximum X value that can be generated
     */
    void _setMaxXValue(const double value);

    /*
     * Set the minimum Y value that can be generated
     */
    void _setMinYValue(const double value);

    /*
     * Set the maximum Y value that can be generated
     */
    void _setMaxYValue(const double value);

    /*
     * Set the distributions to draw numbers from
     */
    void _setDistributions(void);

    /*
     * Generate a random number in [_minX, _maxX)
     */
    double getRandomX(void);

    /*
     * Generate a random number in [_minY, _maxY)
     */
    double getRandomY(void);

    /*
     * Check whether a number is accepted as a sample drawn from a distribution described by func
     */
    bool isAccepted(const double xVal, const double yVal, const double (*func)(double));

  protected:
    /*
     * Minimum X value that can be generated
     */
    double _minX{0};

    /*
     * Maximum X value that can be generated
     */
    double _maxX{0};

    /*
     * Minimum Y value that can be generated
     */
    double _minY{0};

    /*
     * Maximum Y value that can be generated
     */
    double _maxY{0};

  private:
    /*
     * Mersenne Twister engine
     */
    std::mt19937 _gen;

    /*
     * Uniform distribution to draw numbers from when getRandomX() is called
     */
    std::uniform_real_distribution<double> _xDistribution;

    /*
     * Uniform distribution to draw numbers from when getRandomY() is called
     */
    std::uniform_real_distribution<double> _yDistribution;
};

#endif // MCGENERATOR_HPP
