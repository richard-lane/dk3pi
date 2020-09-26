#ifndef MINUIT_FITTER_BASE_H
#define MINUIT_FITTER_BASE_H

#include "BaseFitter.h"
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
     * Perform a fit using Minuit
     *
     * Uses _parameters to minimise _fitFcn; populates this->min, this->fitParams and this->statistic
     *
     * does not set _bestFitFunction (as we don't yet know what to do with our fit
     * parameters); this should be set by child class implementations of this->fit().
     */
    virtual void fit();

    /*
     * Given a vector representing the covariance between a set of parameters, find the covariance matrix using the
     * errors in fitParams.fitParamErrors and convert it to a TMatrixD
     */
    TMatrixD covarianceVector2CorrelationMatrix(const std::vector<double>& covarianceVector);

    /*
     * Save a plot of our data to file
     *
     * Must specify parameters for drawing a legend.
     *
     * Intended to be overridden for plotting best-fit
     */
    virtual void saveFitPlot(const std::string&          plotTitle,
                             const std::string&          path,
                             const util::LegendParams_t* legendParams = nullptr);

    /*
     * Fix parameters by name to their current values.
     */
    void fixParameters(const std::vector<std::string>& fixParams);

    /*
     * Free parameters by name
     */
    void freeParameters(const std::vector<std::string>& freeParams);

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

    /*
     * Pointer to the Minuit FCN used to perform the fit
     */
    std::unique_ptr<MyBaseFcn> _fitFcn = nullptr;

  protected:
    /*
     * Calls parent constructor
     */
    MinuitFitterBase(const FitData_t& fitData);

    /*
     * Helper function to store the result from a Minuit minimisation in this class' fitParams
     */
    void _storeMinuitFitParams();

    /*
     * Set fit parameters and their names
     */
    void _setParams(const std::vector<std::string>& names,
                    const std::vector<double>&      values,
                    const std::vector<double>&      errors);

    /*
     * Parameters that Minuit minimises
     *
     * Should be populated by the parent class.
     */
    std::unique_ptr<ROOT::Minuit2::MnUserParameters> _parameters = nullptr;
};

#endif // MINUIT_FITTER_BASE_H