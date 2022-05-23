#ifndef MULTI_BIN_FIT_H
#define MULTI_BIN_FIT_H

#include "ConstrainedFitter.h"

namespace CharmFitter
{

/*
 * Minuit2-compatible class for fitting charm decay lifetimes
 *
 * Performs fits to the decay times in the four bins in parallel
 *
 * Fits by integrating over time bins
 */
class MultiBinFitFcn : public ROOT::Minuit2::FCNBase
{
  public:
    /*
     * Creates the fit fcns
     */
    MultiBinFitFcn(const std::vector<double>& bin1Data,
                   const std::vector<double>& bin2Data,
                   const std::vector<double>& bin3Data,
                   const std::vector<double>& bin4Data,
                   const std::vector<double>& times,
                   const std::vector<double>& errors,
                   const IntegralOptions_t&   integralOptions);

    /*
     * Creates a model given a vector of {x, y, r1, r2, r3, r4, z_im1, z_im2, z_im3, z_im4, z_re1, z_re2, z_re3, z_re4,
     * width};
     *
     * Returns chi squared value between the model given our params and the measured data, without constraining x and y
     * using the world average data
     *
     */
    virtual double operator()(const std::vector<double>& parameters) const override;

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

  private:
    /*
     * Bin edges + the form of the efficiency function
     */
    IntegralOptions_t _integralOptions;
};

/*
 * Fit ratio of D->K3Pi rates to the equation involving r_D, x, y, etc. across the four bins
 */
class MultiBinFitter
{
  public:
    /*
     * Tell the fitter about the bin limits and what initial guesses of our parameters and errors to make
     */
    MultiBinFitter(const std::vector<double>&    binLimits,
                   const std::array<double, 15>& initialValues,
                   const std::array<double, 15>& initialErrors);

    /*
     * Fit our data to (r^2 + r(yRZ + xImZ)*Gamma*t + (x2 +y2)/4 (Gamma*t)2)
     *
     */
    FitResults_t fit(const std::function<double(double)>& efficiency);
};

} // namespace CharmFitter

#endif // MULTI_BIN_FIT_H
