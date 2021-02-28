#ifndef BASEFITTER_H
#define BASEFITTER_H

#include "FitterUtils.h"
#include "util.h"

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>

namespace CharmFitter
{

/*
 * Base class for fitting a polynomial using Minuit2
 *
 * Does not contain the necessary operator() method needed for Minuit- child classes should define this.
 */
class CharmBaseFcn : public ROOT::Minuit2::FCNBase
{
  public:
    /*
     * Boilerplate that Minuit requires
     */
    virtual double      Up() const { return theErrorDef; };
    std::vector<double> measurements() const { return theMeasurements; };
    std::vector<double> positions() const { return thePositions; };
    std::vector<double> variances() const { return theMVariances; };
    void                setErrorDef(double def) { theErrorDef = def; };

    std::vector<double> theMeasurements;
    std::vector<double> thePositions;
    std::vector<double> theMVariances;
    double              theErrorDef;

    virtual double operator()(const std::vector<double>& parameters) const = 0;

  protected:
    /*
     * Sets the various private members
     */
    CharmBaseFcn(const std::vector<double>& data,
                 const std::vector<double>& times,
                 const std::vector<double>& errors,
                 const IntegralOptions_t&   integralOptions);

    /*
     * Bin edges + the form of the efficiency function
     */
    IntegralOptions_t _integralOptions;
};

/*
 * Virtual base class that other fitters inherit from
 *
 * Defines public APIs for adding data to the histograms + performing the fit, etc.
 *
 */
class CharmFitterBase
{
  public:
    /*
     * Tell the fitter which bin limits to use
     */
    CharmFitterBase(const std::vector<double>& binLimits);

    /*
     * Add a RS event that decayed at the specified time
     *
     * Defined as a method rather than in the constructor as we might have to add points to the histogram one by one if
     * there is too much data to store in an array
     */
    void addRSPoint(const double time, const double weight = 1);

    /*
     * Add a WS event that decayed at the specified time
     *
     * Defined as a method rather than in the constructor as we might have to add points to the histogram one by one if
     * there is too much data to store in an array
     */
    void addWSPoint(const double time, const double weight = 1);

    /*
     * Add a collection of RS events that decayed at the specified times
     *
     * Defined as a method rather than in the constructor as we might have to add points to the histogram one by one if
     * there is too much data to store in an array
     */
    template <typename Container> void addRSPoints(const Container& times, const Container& weights)
    {
        for (size_t i = 0; i < times.size(); ++i) {
            addRSPoint(times[i], weights[i]);
        }
    }

    /*
     * Add a collection of WS events that decayed at the specified times
     *
     * Defined as a method rather than in the constructor as we might have to add points to the histogram one by one if
     * there is too much data to store in an array
     */
    template <typename Container> void addWSPoints(const Container& times, const Container& weights)
    {
        for (size_t i = 0; i < times.size(); ++i) {
            addWSPoint(times[i], weights[i]);
        }
    }

    /*
     * Perform a fit to the data that has been provided
     *
     * Must be implemented in a child class: likely uses _performFit helper method
     */
    virtual FitResults_t fit(const std::function<double(double)>& integralOptions) = 0;

    /*
     * Fix a parameter by name to its current value
     */
    void fixParameter(const std::string& paramName);

    /*
     * Fix parameters by name to their current values.
     */
    template <typename Container> void fixParameters(const Container& params)
    {
        for (const std::string& p : params) {
            fixParameter(p);
        }
    }

    /*
     * Free a parameter by name
     */
    void freeParameter(const std::string& paramName);

    /*
     * Free parameters by name
     */
    template <typename Container> void freeParameters(const Container& params)
    {
        for (const std::string& p : params) {
            freeParameter(p);
        }
    }

    /*
     * Set the value of a parameter
     */
    void setParameter(const std::string& name, const double value) { _parameters->SetValue(name, value); };

    /*
     * Get the centre of each bin
     */
    const std::vector<double>& getBinCentres(void) const { return _binCentres; }

    /*
     * Get the width of each bin
     */
    const std::vector<double>& getBinWidths(void) const { return _binWidths; }

    /*
     * Get the bin content from the RS histogram
     */
    const std::vector<double>& getRSBinContent(void) const { return _rsCounts; }

    /*
     * Get the bin content from the WS histogram
     */
    const std::vector<double>& getWSBinContent() const { return _wsCounts; }

    /*
     * Get the RS/WS ratio in each bin
     */
    std::vector<double> ratios(void) const;

    /*
     * Get the errors on the RS/WS ratio in each bin
     */
    std::vector<double> errors(void) const;

  protected:
    /*
     * Perform a fit
     *
     * Uses _parameters to minimise _fitFcn; populates this->min, this->fitParams and this->statistic
     *
     * does not set _bestFitFunction (as we don't yet know what to do with our fit
     * parameters); this should be set by child class implementations of this->fit().
     *
     * NB: Minuit2 is weird; the data to fit to has to be provided in the Minuit2 FCN constructor, so this must be done
     * in the child class implementation of fit()
     */
    FitResults_t _performFit();

    /*
     * Set fit parameters' values and names
     *
     * Takes vectors rather than arrays to enable us to pass variable numbers of param names/values/errors for the
     * different fit parametrisations
     */
    void _setParams(const std::vector<std::string>& names,
                    const std::vector<double>&      values,
                    const std::vector<double>&      errors);

    /*
     * Given a vector representing the covariance between a set of parameters, find the covariance matrix using the
     * errors in fitParams.fitParamErrors and convert it to a TMatrixD
     */
    TMatrixD _covarianceVector2CorrelationMatrix(const std::vector<double>& covarainceVector,
                                                 const std::vector<double>& fitErrors);

    /*
     * Number of RS points in each bin
     *
     * Double precision as points are likely to be weighted
     */
    std::vector<double> _rsCounts{};

    /*
     * Number of WS points in each bin
     *
     * Double precision as points are likely to be weighted
     */
    std::vector<double> _wsCounts{};

    /*
     * Parameters that Minuit uses for the minimisation
     *
     * Should be populated by the parent class.
     */
    std::unique_ptr<ROOT::Minuit2::MnUserParameters> _parameters = nullptr;

    /*
     * Pointer to the Minuit FCN used to perform the fit
     */
    std::unique_ptr<CharmBaseFcn> _fitFcn = nullptr;

    size_t              _numBins{0};
    std::vector<double> _binEdges{};
    std::vector<double> _binWidths{};
    std::vector<double> _binCentres{};
};

/*
 * Perform a scan of two named parameters over the specified ranges
 *
 * TODO make this return a custom struct which better represents a 2d array
 * Currently returns a vector of shape (range1.size(), range2.size())
 */
std::vector<std::vector<double>> twoDParamScan(CharmFitterBase&                     Fitter,
                                               const std::function<double(double)>& efficiency,
                                               const std::string&                   param1,
                                               const std::string&                   param2,
                                               const std::vector<double>&           range1,
                                               const std::vector<double>&           range2);

/*
 * Save a plot of the data stored in this fitter, with the given best fit function superimposed
 *
 * Plots points at the centre of each bin, with error bars covering the bin width
 */
void savePlot(const CharmFitterBase&      Fitter,
              TF1&                        bestFitFunction,
              const std::string&          path,
              const std::string&          title,
              const util::LegendParams_t& legend);

} // namespace CharmFitter

#endif // BASEFITTER_H
