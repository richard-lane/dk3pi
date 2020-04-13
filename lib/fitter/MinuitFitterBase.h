#ifndef MINUIT_FITTER_BASE_H
#define MINUIT_FITTER_BASE_H

#include "BaseFitter.h"
#include "MinuitFitter.h"
#include "util.h"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "TF1.h"

/*
 * Base class for fitters that use Minuit
 */
class MinuitFitterBase : public BaseFitter
{
  public:
    /*
     * Perform a fit using Minuit, possibly holding some parameters fixed as specified in fixParams.
     *
     * Uses _parameters to minimise _fitFcn; populates this->min, this->fitParams and this->statistic
     *
     * does not set _bestFitFunction (as we don't yet know what to do with our fit
     * parameters); this should be set by child class implementations of this->fit().
     *
     */
    virtual void fit(const std::vector<size_t>& fixParams);

    /*
     * Set our minimiser using this class' _fitFcn and _parameters
     */
    void setMigrad();

    /*
     * Reset our minimiser
     */
    void resetMigrad();

    /*
     * Set fit parameters and their names
     */
    void setParams(const std::vector<std::string>& names,
                   const std::vector<double>&      values,
                   const std::vector<double>&      errors);

    /*
     * Given a vector representing the covariance between a set of parameters, find the covariance matrix using the
     * errors in fitParams.fitParamErrors and convert it to a TMatrixD
     *
     * Must tell this function which parameters were fixed; {0, 2, 3} means params 0, 2 and 3 were fixed.
     */
    TMatrixD covarianceVector2CorrelationMatrix(const std::vector<double>& covarianceVector,
                                                const std::vector<size_t>& fixedParams);

    /*
     * Save a plot of our data and fit to file
     *
     * Must specify parameters for drawing a legend.
     */
    void saveFitPlot(const std::string&          plotTitle,
                     const std::string&          path,
                     const util::LegendParams_t* legendParams = nullptr);

    /*
     * Our best-fit function
     *
     * Is not set by this class' fit function (as we don't yet know what to do with our fit parameters); should be set
     * by child class implementations of this->fit().
     */
    std::unique_ptr<TF1> bestFitFunction = nullptr;

    /*
     * fit status etc
     */
    std::unique_ptr<ROOT::Minuit2::FunctionMinimum> min = nullptr;

  protected:
    /*
     * Calls parent constructor
     */
    MinuitFitterBase(const FitData_t& fitData);

    /*
     * Helper function to store the result from a Minuit minimisation in this class' fitParams
     */
    void _storeMinuitFitParams(const std::vector<size_t>& fixIndices);

    /*
     * Pointer to the Minuit FCN used to perform the fit
     */
    std::unique_ptr<BasePolynomialFcn> _fitFcn = nullptr;

    /*
     * Parameters that Minuit minimises
     *
     * Should be populated by the parent class.
     */
    std::unique_ptr<ROOT::Minuit2::MnUserParameters> _parameters = nullptr;

    /*
     * Minuit2 minimiser
     */
    std::unique_ptr<ROOT::Minuit2::MnMigrad> _migrad = nullptr;
};

#endif // MINUIT_FITTER_BASE_H
