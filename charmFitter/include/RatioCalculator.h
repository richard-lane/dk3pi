/*
 * Bin decay times into specified time bins and calculate the ratio of CF to DCS decays in each.
 */
#ifndef RATIOCALCULATOR_H
#define RATIOCALCULATOR_H

#include <utility>
#include <vector>

class RatioCalculator
{
  public:
    /*
     * Constructor
     *
     * Check that bin limits are sorted and populate decay time and bin data.
     */
    RatioCalculator(const std::vector<size_t> &denominator,
                    const std::vector<size_t> &numerator,
                    const std::vector<double> &binLimits);

    /*
     * Bin the times passed to this class and find their ratio.
     * Populates this class' ratio and error attributes.
     */
    void calculateRatios(void);

    /*
     * Ratio of the number of DCS to CF decays in each bin
     */
    std::vector<double> ratio{};

    /*
     * Error in the ratio of number of DCS to CF decays
     */
    std::vector<double> error{};

    /*
     * Centres and widths of the time bins
     */
    std::vector<double> binCentres{};
    std::vector<double> binWidths{};

  private:
    std::vector<size_t> _denominator{};
    std::vector<size_t> _numerator{};

    /*
     * A vector containing every bin edge, from the left edge of the lowest bin to the right edge of the highest
     */
    std::vector<double> _binLimits{};

    /*
     * Number of time bins
     */
    size_t _numBins{0};
};

#endif // RATIOCALCULATOR_H
