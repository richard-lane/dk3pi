/*
 * Fit binned data
 */
#ifndef FITTER_H
#define FITTER_H

#include <memory>
#include <vector>

#include "DecaySimulator.h"

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
     * Fit our data to a second order polynomial, using ROOT's builtin "pol2" function
     * Sets fitParams attribute and allocates memory to _plot.
     *
     * Fit options can be passed via the options argument.
     *
     * Populates fitParams
     */
    void fitUsingRootBuiltinPol2(const std::string& options = "");

    /*
     * Fit our data to the equation we expect to see, using ROOT's builtin TGraph fitter.
     * At the moment this is just a second order polynomial.
     *
     * minTime and maxTime define the range over which the function is defined.
     *
     * Populates fitParams
     */
    void fitUsingRootCustomFcn(const double minTime, const double maxTime, const std::string& options = "");

    /*
     * Fit our data to a second-order polynomial a + bt + ct^2 using Minuit2 and the chi-squared method.
     *
     * The user should provide an initial guess at the parameters and their errors
     *
     * Populates fitParams
     */
    void fitUsingMinuit2ChiSq(const std::vector<double>& initialParams, const std::vector<double>& initialErrors);

    /*
     * Given a vector representing the covariance between a set of parameters,  find the covariance matrix using the
     * errors in fitParams.fitParamErrors and convert it to a TMatrixD
     */
    TMatrixD covarianceVector2CorrelationMatrix(const std::vector<double>& covarianceVector);

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
