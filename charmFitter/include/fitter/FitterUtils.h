/*
 * Utility structs and stuff for fitting
 * 
 * Desperately needs renaming
 */
#ifndef FITTERUTILS_H
#define FITTERUTILS_H

#include <vector>

#include "Minuit2/FCNBase.h"
#include "TMatrixD.h"

/*
 * The data needed to perform the fit
 * Bin limits, our data and its error
 *
 * Parameters are all initialised to zero by default and should be set by the user.
 *
 * Bin limits should be sorted
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
    FitData(const std::vector<double>& myBinLimits,
            const std::vector<double>& myData,
            const std::vector<double>& myErrors);

    std::vector<double> binLimits;
    std::vector<double> data;
    std::vector<double> errors;
    size_t              numBins{0};

    /*
     * Centre of each bin
     */
    std::vector<double> binCentres;

} FitData_t;

/*
 * The parameters describing how a function fits data
 */
typedef struct FitResults {

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
    std::unique_ptr<TMatrixD> correlationMatrix{};

    /*
     * Fit statistic, e.g. chi squared or likelihood
     */
    double fitStatistic{0.0};
} FitResults_t;

/*
 * Struct encapsulating the additional data a fitter needs to know about if it is to integrate over the bins in the chi
 * squared model.
 *
 * If integrate=false, the other options are moot
 */
typedef struct IntegralFitOptions {
    // Initialises everything to 0
    IntegralFitOptions(){};

    // Initialises to the provided values
    IntegralFitOptions(const bool                integrate,
                       const double              width,
                       const std::vector<double> binLimits,
                       const double              efficiencyTimescale)
        : binLimits(binLimits), width(width), efficiencyTimescale(efficiencyTimescale), integrate(integrate){};

    std::vector<double> binLimits{0};           // Our data's bin limits
    double              width{0};               // Decay width
    double              efficiencyTimescale{0}; // Characteristic time parametrising the efficiency function

    bool integrate{false};
} IntegralOptions_t;

/*
 * Base class for fitting a polynomial using Minuit2
 *
 * Does not contain the necessary operator() method needed for Minuit- child classes should define this.
 */
class MyBaseFcn : public ROOT::Minuit2::FCNBase
{
  public:
    /*
     * Maybe we need a destructor
     */
    ~MyBaseFcn();

    /*
     * Returns theErrorDef
     */
    virtual double Up() const;

    /*
     * Return the values of the function
     */
    std::vector<double> measurements() const;

    /*
     * Return the positions where our datapoints are evaluated
     */
    std::vector<double> positions() const;

    /*
     * Return variances
     */
    std::vector<double> variances() const;

    /*
     * Sets whatever theErrorDef is
     */
    void setErrorDef(double def);

    /*
     * Measurements; our datapoints
     */
    std::vector<double> theMeasurements;

    /*
     * Where our datapoints are evaluated
     */
    std::vector<double> thePositions;

    /*
     * Errors
     */
    std::vector<double> theMVariances;

    /*
     * Used to define errors. Used internally by Minuit so just don't worry about it
     */
    double theErrorDef;

  protected:
    /*
     * Sets the various private members
     */
    MyBaseFcn(const std::vector<double>& data,
              const std::vector<double>& times,
              const std::vector<double>& errors,
              const IntegralOptions_t&   integralOptions);

    /*
     * Whether to integrate + options if so
     */
    IntegralOptions_t _integralOptions;
};

#endif // FITTERUTILS_H
