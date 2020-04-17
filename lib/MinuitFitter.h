/*
 * Fit binned data with Minuit
 */
#ifndef MINUIT_FITTER_H
#define MINUIT_FITTER_H

#include <memory>
#include <vector>

#include "DecaySimulator.h"

#include "Minuit2/FCNBase.h"
#include "TGraphErrors.h"
#include "TMatrixD.h"

/*
 * Base class for fitting a polynomial using Minuit2
 *
 * Does not contain the necessary operator() method needed for Minuit- child classes should define this.
 */
class BasePolynomialFcn : public ROOT::Minuit2::FCNBase
{
  public:
    /*
     * Sets the various private members
     */
    BasePolynomialFcn(const std::vector<double>& data,
                      const std::vector<double>& times,
                      const std::vector<double>& errors);

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
                       const std::vector<double>& binLimits,
                       const double               width);

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
     * Bin limits that were used in performing this fit
     */
    const std::vector<double> _binLimits{};

    /*
     * Decay width
     */
    const double _width{};
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

#endif // MINUIT_FITTER_H
