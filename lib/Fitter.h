/*
 * Fit binned data
 */
#ifndef FITTER_H
#define FITTER_H

#include <memory>
#include <vector>

#include "DecaySimulator.h"
#include "MinuitFitter.h"
#include "util.h"

#include "Minuit2/FunctionMinimum.h"
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
 * Enum telling the fitter what algorithm to use, if using Minuit2
 */
typedef enum FitAlgorithm { ChiSquared } FitAlgorithm_t;

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
     * Fit our data to the equation we expect to see, using ROOT's builtin TGraph fitter.
     * At the moment this is just a second order polynomial.
     *
     * Allocates memory to _plot
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
     * Parameters are {x, y, r, z_im, z_re, width}
     *
     * FitMethod tells the fitter whether to use chi squared or maximum likelihood (max likelihood isn't actually
     * implemented)
     *
     * Allocates memory to _plot and _bestFitPlot
     *
     * Populates fitParams
     */
    void fitUsingMinuit(const std::vector<double>& initialParams,
                        const std::vector<double>& initialErrors,
                        const FitAlgorithm_t&      FitMethod);

    /*
     * Fit our data to a second-order polynomial r2 + r(yRZ + xImZ)Gt + (x2+y2)(Gt)2/4 using Minuit2 and the chi-squared
     * method.
     *
     * The user should provide an initial guess at the parameters and their errors
     *
     * Allocates memory to _plot and _bestFitPlot
     *
     * Populates fitParams
     */
    void detailedFitUsingMinuit(const std::vector<double>& initialParams,
                                const std::vector<double>& initialErrors,
                                const FitAlgorithm_t&      FitMethod);

    /*
     * Given a vector representing the covariance between a set of parameters,  find the covariance matrix using the
     * errors in fitParams.fitParamErrors and convert it to a TMatrixD
     */
    TMatrixD covarianceVector2CorrelationMatrix(const std::vector<double>& covarianceVector);

    /*
     * Scan the ith parameter as defined in fitParams.fitParams
     *
     * Cannot have more than 100 points due to some of Minuit's internal limitations.
     * By default runs a scan from +-2sigma, but can optionally be specified by setting low and high.
     *
     * Populates parameterScan
     */
    void chiSqParameterScan(const size_t i, const size_t numPoints, const double low = 0., const double high = 0.);

    /*
     * Scan the i and jth parameters between the specified limits
     *
     * Populates twoDParameterScan
     */
    void twoDParamScan(const size_t i,
                       const size_t j,
                       const size_t iPoints,
                       const size_t jPoints,
                       const double iLow,
                       const double iHigh,
                       const double jLow,
                       const double jHigh);

    /*
     * Save a plot of our data and fit to file
     *
     * If plotting from a minuit fit, must specify parameters for drawing a legend.
     */
    void saveFitPlot(const std::string&          plotTitle,
                     const std::string&          path,
                     const util::LegendParams_t* legendParams = nullptr);

    /*
     * Parameters describing the fit
     */
    FitResults_t fitParams;

    /*
     * Vector of pairs describing a parameter scan
     */
    std::vector<std::pair<double, double>> parameterScan;

    /*
     * Vector of tuples describing a 2d parameter scan
     *
     * Scans parameters i and j to find chi squared values; result is a vector of (i_value, j_value, chi_squared)
     */
    std::vector<std::vector<double>> twoDParameterScan;

    /*
     * ROOT TGraph object used for holding input data. Can also be used to perform one of ROOT's builtin fits.
     */
    std::unique_ptr<TGraphErrors> plot = nullptr;

    /*
     * ROOT TGraph object that used for representing the best-fit of our data, in the case that ROOT is not used to
     * perform a builtin fit (e.g. when Minuit2 is used directly.)
     */
    std::unique_ptr<TGraph> bestFitPlot = nullptr;

    /*
     * The test statistic that is optimised when running our test, e.g. chi squared or a likelihood
     */
    std::unique_ptr<double> statistic = nullptr;

    /*
     * fit status etc
     */
    std::unique_ptr<ROOT::Minuit2::FunctionMinimum> min= nullptr;

  private:
    /*
     * Helper function to store the attributes from a Minuit2 FunctionMinimum in this class' fitParams
     */
    void _storeMinuitFitParams(const ROOT::Minuit2::FunctionMinimum& min);

    /*
     * The data used to make the fit
     */
    FitData_t _fitData;

    /*
     * Pointer to the Minuit FCN used to perform the fit
     */
    std::unique_ptr<BasePolynomialFcn> _fitFcn = nullptr;
};

#endif // FITTER_H
