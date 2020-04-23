/*
 * Fit data with Minuit
 *
 * Minuit uses a special class to perform the fit; the ones I might want to use are defined here.
 */
#ifndef MINUIT_FITTER_H
#define MINUIT_FITTER_H

#define WORLD_AVERAGE_X (0.0039183)
#define WORLD_AVERAGE_X_ERR (0.0011489)
#define WORLD_AVERAGE_Y (0.0065139)
#define WORLD_AVERAGE_Y_ERR (0.00064945)
#define X_Y_CORRELATION (-0.301)

#include <memory>
#include <vector>

#include "DecaySimulator.h"

#include "Minuit2/FCNBase.h"
#include "TGraphErrors.h"
#include "TMatrixD.h"

/*
 * Struct encapsulating the additional data a fitter needs to know about if it is to integrate over the bins in the chi
 * squared model.
 */
typedef struct IntegralFitOptions {
    // Initialises everything to 0
    IntegralFitOptions(){};

    // Initialises to the provided values
    IntegralFitOptions(const double              width,
                       const std::vector<double> binLimits,
                       const double              tolerance      = 1e-12,
                       const size_t              maxRefinements = 10)
        : tolerance(tolerance), maxRefinements(maxRefinements), binLimits(binLimits), width(width){};

    double tolerance{0};      // Tolerance to use for numerical integration
    size_t maxRefinements{0}; // Maximum number of refinements to make

    std::vector<double> binLimits{0}; // Our data's bin limits
    double              width{0};     // Decay width
} IntegralOptions_t;

/*
 * Base class for fitting a polynomial using Minuit2
 *
 * Does not contain the necessary operator() method needed for Minuit- child classes should define this.
 */
class BasePolynomialFcn : public ROOT::Minuit2::FCNBase
{
  public:
    /*
     * Maybe we need a destructor
     */
    ~BasePolynomialFcn();

    /*
     * Returns _theErrorDef (i still don't really know what this is )
     */
    virtual double Up() const;

    /*
     * Return the values of the function
     */
    std::vector<double> measurements() const;

    /*
     * Return the positions where our datapoints are evaluated...?
     */
    std::vector<double> positions() const;

    /*
     * Return variances
     */
    std::vector<double> variances() const;

    /*
     * Sets whatever _theErrorDef is
     */
    void setErrorDef(double def);

    /*
     * Measurements; our datapoints
     */
    std::vector<double> theMeasurements;

    /*
     * Where our datapoints are evaluated ...?
     */
    std::vector<double> thePositions;

    /*
     * Errors
     */
    std::vector<double> theMVariances;

    /*
     * I don't know
     */
    double theErrorDef;

  protected:
    /*
     * Sets the various private members
     */
    BasePolynomialFcn(const std::vector<double>& data,
                      const std::vector<double>& times,
                      const std::vector<double>& errors);
};

/*
 * Class fitting to a polynomial using chi squared
 *
 * Fits by integrating over the bins
 */
class PolynomialChiSqFcn : public BasePolynomialFcn
{
  public:
    /*
     * Calls parent constructor
     */
    PolynomialChiSqFcn(const std::vector<double>& data,
                       const std::vector<double>& times,
                       const std::vector<double>& errors,
                       const IntegralOptions_t&   integralOptions);

    ~PolynomialChiSqFcn();

    /*
     * Creates a model given a vector of three parameters.
     *
     * Returns the chi squared value between the model given our parameters and the measured data.
     *
     */
    virtual double operator()(const std::vector<double>& parameters) const;

  private:
    /*
     * options to use while integrating
     */
    IntegralOptions_t _integralOptions;
};

/*
 * Class fitting to a polynomial using chi squared
 *
 * Does not integrate over the bins
 */
class PolynomialChiSqFcnNoIntegral : public BasePolynomialFcn
{
  public:
    /*
     * Calls parent constructor
     */
    PolynomialChiSqFcnNoIntegral(const std::vector<double>& data,
                                 const std::vector<double>& times,
                                 const std::vector<double>& errors);

    ~PolynomialChiSqFcnNoIntegral();

    /*
     * Creates a model given a vector of three parameters.
     *
     * Returns the chi squared value between the model given our parameters and the measured data.
     *
     */
    virtual double operator()(const std::vector<double>& parameters) const;
};

/*
 * Class fitting to a detailed polynomial (r^2 + r(yRZ + xImZ)*Gamma*t + (x2 +y2)/4 (Gamma*t)2)
 */
class DetailedPolynomialChiSqFcn : public BasePolynomialFcn
{
  public:
    /*
     * Calls parent constructor
     */
    DetailedPolynomialChiSqFcn(const std::vector<double>& data,
                               const std::vector<double>& times,
                               const std::vector<double>& errors);

    ~DetailedPolynomialChiSqFcn();

    /*
     * Creates a model given a vector of {x, y, r, z_im, z_re, width};
     *
     * Returns chi squared value between the model given our params and the measured data
     */
    virtual double operator()(const std::vector<double>& parameters) const;
};

/*
 * Class fitting to a detailed polynomial (r^2 + r(yRZ + xImZ)*Gamma*t + (x2 +y2)/4 (Gamma*t)2), with constraint on X
 * and Y
 */
class DetailedChiSqConstrainXYFcn : public BasePolynomialFcn
{
  public:
    /*
     * Calls parent constructor
     */
    DetailedChiSqConstrainXYFcn(const std::vector<double>& data,
                                const std::vector<double>& times,
                                const std::vector<double>& errors);

    ~DetailedChiSqConstrainXYFcn();

    /*
     * Creates a model given a vector of {x, y, r, z_im, z_re, width};
     *
     * Returns chi squared value between the model given our params and the measured data, constraining X and Y to the
     * values
     */
    virtual double operator()(const std::vector<double>& parameters) const;
};

#endif // MINUIT_FITTER_H
