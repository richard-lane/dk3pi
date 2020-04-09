#ifndef MINUIT_FITTER_BASE_H
#define MINUIT_FITTER_BASE_H

#include "BaseFitter.h"
#include "MinuitFitter.h"
#include "util.h"

#include "Minuit2/FunctionMinimum.h"

/*
 * Base class for fitters that use Minuit
 */
class MinuitFitter : public BaseFitter
{
  public:
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
     * ROOT TGraph object that used for representing the best-fit of our data
     * Maybe this should really be a TF1 TODO
     */
    std::unique_ptr<TGraph> bestFitPlot = nullptr;

    /*
     * fit status etc
     */
    std::unique_ptr<ROOT::Minuit2::FunctionMinimum> min = nullptr;

  protected:
    /*
     * Calls parent constructor
     */
    MinuitFitter(const FitData_t& fitData);

    /*
     * Helper function to store the attributes from a Minuit2 FunctionMinimum in this class' fitParams
     */
    void _storeMinuitFitParams(const ROOT::Minuit2::FunctionMinimum& min);

    /*
     * Pointer to the Minuit FCN used to perform the fit
     */
    std::unique_ptr<BasePolynomialFcn> _fitFcn = nullptr;
};

#endif // MINUIT_FITTER_BASE_H