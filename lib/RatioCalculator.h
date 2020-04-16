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
     * Check that bin limits are sorted and populate decay time and bin data.
     */
    RatioCalculator(const std::vector<size_t> &cfDecayCounts,
                    const std::vector<size_t> &dcsDecayCounts,
                    const std::vector<double> &binLimits);

    /*
     * Find the ratio of each element of two vectors, and its associated error assuming both vectors obey Poisson
     * statistics.
     * Returns a vector of std::pairs (ratio, error).
     *
     * Public to allow for unit testing
     */
    std::vector<std::pair<double, double>> findRatioAndError(const std::vector<size_t> &numerator,
                                                             const std::vector<size_t> &denominator);

    /*
     * Bin the times passed to this class and find their ratio.
     * Populates this class' ratio and error attributes.
     *
     * This is what a user should call to populate this class with appropriate data.
     */
    void calculateRatios(void);

    /*
     * Save the numbers of DCS and CF points in each bin to file
     */
    void findNumPointsPerBin(const std::string &path);

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
    /*
     * Vectors of CF and DCS decay times
     */
    std::vector<size_t> _cfDecayCounts{};
    std::vector<size_t> _dcsDecayCounts{};

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
