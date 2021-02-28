/*
 * Utility structs and stuff for fitting
 *
 * Desperately needs renaming
 */
#ifndef FITTERUTILS_H
#define FITTERUTILS_H

#include <vector>
#include <memory>

#include <Minuit2/FCNBase.h>
#include <TF1.h>
#include <TMatrixD.h>

namespace CharmFitter
{

/*
 * The parameters describing how a function fits data
 */
typedef struct FitResults {
    /*
     * Fit statistic, e.g. chi squared or likelihood
     */
    double fitStatistic{0.0};

    /*
     * Vector of the fit parameters
     */
    std::vector<double> fitParams{};

    /*
     * Vector of the errors in our fit parameters
     */
    std::vector<double> fitParamErrors{};

    /*
     * Pointer to correlation matrix of fit parameters
     */
    TMatrixD correlationMatrix{};

    /*
     * Pointer to ROOT best-fit function
     */
    TF1 bestFitFunction{};
} FitResults_t;

/*
 * Struct encapsulating the additional data a fitter needs to know about if it is to integrate over the bins in the chi
 * squared model.
 *
 */
typedef struct IntegralFitOptions {
    /*
     * The bin limits used when binning the datasets into histograms
     */
    const std::vector<double>& binLimits;

    /*
     * The decay-time efficiency function
     */
    const std::function<double(double)>& efficiency;
} IntegralOptions_t;

} // namespace CharmFitter

#endif // FITTERUTILS_H
