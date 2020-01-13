#ifndef DATA_SETS_RATIO_CPP
#define DATA_SETS_RATIO_CPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TROOT.h"

#include "../include/DataSetsRatio.h"

DataSetsRatio::DataSetsRatio(std::vector<size_t> &myNumeratorData,
                             std::vector<size_t> &myDenominatorData,
                             std::vector<double> &myBinCentres,
                             std::vector<double> &myBinErrors)
{
    verifyInputs(myNumeratorData, myDenominatorData, myBinCentres, myBinErrors);

    // Set bin limits
    binCentres = myBinCentres;
    binErrors  = myBinErrors;

    numeratorData   = myNumeratorData;
    denominatorData = myDenominatorData;

    numBins = binCentres.size();

    _setRatios();
}

/*
 * Check that our datasets have the right number of points.
 */
void DataSetsRatio::verifyInputs(std::vector<size_t> &myNumeratorData,
                                 std::vector<size_t> &myDenominatorData,
                                 std::vector<double> &myBinCentres,
                                 std::vector<double> &myBinErrors)
{
    numBins = myBinCentres.size();
    if ((myNumeratorData.size() != numBins || myDenominatorData.size() != numBins) || myBinErrors.size() != numBins) {
        std::cerr << "Incompatible bins and datasets provided; must have the same number of bins and datapoints"
                  << std::endl;
        throw;
    }
}

/*
 * Remove NaN, inf and 0 from our vectors of data.
 */
void DataSetsRatio::_pruneBadRatios()
{
    // Prune NaN and Inf from ratios.
    // Also remove the corresponding error and bin
    numPoints = numBins;
    auto it   = ratios.begin();
    while (it != ratios.end()) {
        // Use TMath::Finite instead of std::is_finite() because for some reason including the ROOT headers makes
        // std::is_finite(-inf) sometimes return true...
        if (!TMath::Finite(*it) || *it == 0.0) {
            size_t index = it - ratios.begin();

            ratios.erase(ratios.begin() + index);
            ratioErrors.erase(ratioErrors.begin() + index);
            binCentres.erase(binCentres.begin() + index);
            binErrors.erase(binErrors.begin() + index);

            numPoints--;

        } else {
            ++it;
        }
    }
}

/*
 * Set the ratio of our numerator and denominator's points in each bin
 */
void DataSetsRatio::_setRatios()
{
    // Reassign our vector of ratios to -1
    ratios.assign(numBins, -1);
    ratios.shrink_to_fit();

    // Unintelligently divide our elements
    // A good implementation would use std::transform but this is fine
    for (size_t i = 0; i < numBins; ++i) {
        ratios[i] = (double)numeratorData[i] / (double)denominatorData[i];
    }

    // Find also the errors in our ratios.
    _setRatioErrors();

    _pruneBadRatios();
}

/*
 * Assuming the error in our counts is sqrt(count) and that the fractional error in a ratio is given by adding count
 * errors in quadrature, find the error in a ratio.
 *
 * Sets errors for ratio=0 to 0
 */
double DataSetsRatio::ratioError(const double &ratio, const size_t &numeratorCounts, const size_t &denominatorCounts)
{
    // If our ratio is zero, our error should also be
    if (ratio == 0) {
        return 0;
    }

    // Otherwise our ratio is found by assuming delta(Counts) = sqrt(Counts)
    return std::sqrt(((double)numeratorCounts + (double)denominatorCounts) /
                     ((double)numeratorCounts * (double)denominatorCounts)) *
           ratio;
}

/*
 * Set the values of ratioErrors
 */
void DataSetsRatio::_setRatioErrors()
{
    ratioErrors.assign(numBins, 0);
    ratioErrors.shrink_to_fit();

    for (size_t i = 0; i < numBins; ++i) {
        ratioErrors[i] = ratioError(ratios[i], numeratorData[i], denominatorData[i]);
    }
}

/*
 * Plot the ratios of numerator to denominator points in each bin
 * Doesn't actually draw the plot but assigns _ratioPlot to a new TGraphErrors object
 */
void DataSetsRatio::plotBinRatios()
{
    _ratioPlot = new TGraphErrors(numPoints, binCentres.data(), ratios.data(), binErrors.data(), ratioErrors.data());
}

/*
 * Fit a second order polynomial to a plot of ratios in each bin
 *
 * @param draw: whether to draw the graph or just fit to its data.
 */
void DataSetsRatio::fitToData(bool draw, std::string plotTitle)
{
    // Create a plot for our data, giving it a custom global title and titling the x and y axes with time and ratio
    plotBinRatios();
    _ratioPlot->SetTitle((plotTitle + ";time/ns;CF/DCS ratio").c_str());

    if (draw) {
        TCanvas *c = new TCanvas();
        _ratioPlot->Draw("*ap");
    }

    // ROOT is terrible and will often fail to fit when the TGraph contains x-errors, unless the initial fit parameters
    // are reasonably close to the true values. We can get around this by perfoming an initial fit ignoring error bars
    // ("W"), then fitting again.
    TF1 *graph_fit = ((TF1 *)(gROOT->GetFunction("pol2")));

    _ratioPlot->Fit(graph_fit, "WQRN");
    _ratioPlot->Fit(graph_fit);

    // Print a newline after the fitter output to make things easier to read.
    std::cout << std::endl;
}

#endif // DATA_SETS_RATIO_CPP
