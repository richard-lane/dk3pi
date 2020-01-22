/*
 * Monte-Carlo generation of points in a distribution
 */

#ifndef MCGENERATOR_HPP
#define MCGENERATOR_HPP

#include <random>
#include <utility>

/*
 * Base class for generation of points using Monte Carlo
 *
 * Contains methods for generating from a distribution and running accept/reject
 * Intented to be a parent class to useful implementations
 */
class MCGenerator
{
  public:
    MCGenerator(const std::pair<double, double> &xRange = std::make_pair(0.0, 1.0),
                const std::pair<double, double> &yRange = std::make_pair(0.0, 1.0));

    /*
     * Set the range of X values that can be generated
     */
    void setXRange(const std::pair<double, double> &xRange);

    /*
     * Set the range of Y values that can be generated
     */
    void setYRange(const std::pair<double, double> &yRange);

    /*
     * Set the distributions to draw numbers from
     */
    void setDistributions(void);

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
