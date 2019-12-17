#ifndef DATA_SETS_RATIO_HPP
#define DATA_SETS_RATIO_HPP

#include <vector>

#include "TGraph.h"

/*
 * Class for calculating and representing the ratio of two binned datasets
 *
 * A ratio of the two datasets will be calculated, so one is referred to as "numerator" and one as "denominator"
 *
 */
class DataSetsRatio
{
  public:
    explicit DataSetsRatio(std::vector<size_t> &myNumeratorData,
                           std::vector<size_t> &myDenominatorData,
                           std::vector<double> &myBinCentres,
                           std::vector<double> &myBinErrors);

    void _setRatios();
    void fitToData(bool draw, std::string plotTitle);

  private:
    void verifyInputs(std::vector<size_t> &myNumeratorData,
                      std::vector<size_t> &myDenominatorData,
                      std::vector<double> &myBinCentres,
                      std::vector<double> &myBinErrors);
    void _setRatioErrors();
    void _pruneBadRatios();
    void plotBinRatios();

    double ratioError(const double &ratio, const size_t &numeratorCounts, const size_t &denominatorCounts);

    TGraph *_ratioPlot = nullptr;

    // Once we remove any bins that have inf/NaN/zeroes in their data, our number of data points will be smaller than
    // the number of bins.
    size_t numBins{0};
    size_t numPoints{0};

    std::vector<size_t> numeratorData{};
    std::vector<size_t> denominatorData{};

    std::vector<double> ratios{};
    std::vector<double> ratioErrors{};

    std::vector<double> binCentres{};
    std::vector<double> binErrors{};
};

#endif // DATA_SETS_RATIO_HPP