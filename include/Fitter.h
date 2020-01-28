/*
 * Fit binned data
 */
#ifndef FITTER_H
#define FITTER_H

#include <memory>
#include <vector>

#include "TGraphErrors.h"
#include "TMatrixD.h"

/*
 * The data needed to perform the fit
 * Bin centres, bin widths, our data and its error
 *
 * Parameters are all initialised to zero by default and should be set by the user.
 */
typedef struct FitData {

    /*
     * Default constructor sets all attributes to zero
     */
    FitData(void);

    /*
     * Check that all parameters are vectors of the same length, our bins don't overlap and data is all finite.
     * Note that this will not accept a dataset which contains zeros.
     * Sets numPoints
     */
    FitData(const std::vector<double>& myBinCentres,
            const std::vector<double>& myBinWidths,
            const std::vector<double>& myData,
            const std::vector<double>& myErrors);

    std::vector<double> binCentres{0.0};

    /*
     * This is not the total width of the bins; it is the error in the bin, or half the width
     */
    std::vector<double> binErrors{0.0};

    std::vector<double> data{0.0};
    std::vector<double> errors{0.0};

    /*
     * The number of data points that our dataset represents
     * Set by FitData constructor.
     */
    size_t numPoints{0};
} FitData_t;

/*
 * The parameters describing how a polynomial fits data
 */
typedef struct PolynomialFitResults {

    /*
     * Vector of the fit parameters (a0, a1, a2, ...) for a polynomial a0 + a1 * x + a2 * x^2 + ...
     */
    std::vector<double> fitParams{};

    /*
     * Vector of the errors in our fit parameters
     */
    std::vector<double> fitParamErrors{};

    /*
     * Pointer to correlation matrix of fit parameters
     */
    std::unique_ptr<TMatrixD> correlationMatrix{};
} FitResults_t;

/*
 * Class for fitting data to a second-order polynomial
 */
class Fitter
{
  public:
    /*
     * Tell the Fitter the data to be fit.
     */
    Fitter(const FitData_t& fitData);

    /*
     * Fit our data to a second order polynomial.
     * Sets fitParams attribute and allocates memory to _plot.
     *
     * Fit options can be passed via the options argument.
     */
    void pol2fit(const std::string& options = "");

    /*
     * Save a plot of our data and fit to file
     */
    void saveFitPlot(const std::string& plotTitle, const std::string& path);

    /*
     * Parameters describing the fit
     */
    FitResults_t fitParams;

  private:
    /*
     * The data used to make the fit
     */
    FitData_t _fitData;

    /*
     * ROOT TGraph object used for performing the fit
     */
    std::unique_ptr<TGraphErrors> _plot = nullptr;
};

#endif // FITTER_H
