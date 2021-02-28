#ifndef CHARM_UNCONSTRAINED_FITTER_H
#define CHARM_UNCONSTRAINED_FITTER_H

#include "CharmFitterBase.h"

namespace CharmFitter
{

/*
 * Minuit2-compatible class for fitting charm decay lifetimes without constraining
 * mixing parameters x and y
 *
 * Fits by integrating over the bins
 */
class UnconstrainedChiSqFcn : public CharmBaseFcn
{
  public:
    /*
     * Calls parent constructor
     */
    UnconstrainedChiSqFcn(const std::vector<double>& data,
                          const std::vector<double>& times,
                          const std::vector<double>& errors,
                          const IntegralOptions_t&   integralOptions);

    /*
     * Creates a model given a vector of three parameters.
     *
     * Returns the chi squared value between the model given our parameters and the measured data.
     *
     */
    virtual double operator()(const std::vector<double>& parameters) const override;
};

/*
 * Fit ratio of D->K3Pi rates to the equation involving r_D, x, y, etc.; without constraining x and y to their world
 * average values
 */
class UnconstrainedFitter : public CharmFitterBase
{
  public:
    /*
     * Tell the fitter about the raw data and what initial guesses of our parameters and errors to make
     *
     * Parameter ordering is defined by _paramNames
     */
    UnconstrainedFitter(const std::vector<double>&   binLimits,
                        const std::array<double, 6>& initialValues,
                        const std::array<double, 6>& initialErrors);

    /*
     * Fit our data to (r^2 + r(yRZ + xImZ)*Gamma*t + (x2 +y2)/4 (Gamma*t)2)
     *
     */
    FitResults_t fit(const std::function<double(double)>& efficiency) override;

  private:
    /*
     * Parameters used in the fit
     *
     * Useful for indexing + consistency; we can use this to e.g. fix a parameter by name
     */
    const static std::array<std::string, 6> _paramNames;
};

} // namespace CharmFitter

#endif // CHARM_UNCONSTRAINED_FITTER_H
