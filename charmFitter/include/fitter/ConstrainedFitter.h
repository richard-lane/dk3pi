#ifndef CHARM_CONSTRAINED_FITTER_H
#define CHARM_CONSTRAINED_FITTER_H

#include "CharmFitterBase.h"

namespace CharmFitter
{
constexpr double WORLD_AVERAGE_X     = 0.0039183;
constexpr double WORLD_AVERAGE_X_ERR = 0.0011489;
constexpr double WORLD_AVERAGE_Y     = 0.0065139;
constexpr double WORLD_AVERAGE_Y_ERR = 0.00064945;
constexpr double X_Y_CORRELATION     = -0.301;

/*
 * Minuit2-compatible class for fitting charm decay lifetimes, constraining
 * mixing parameters x and y to their world averages
 *
 * Fits by integrating over the bins
 */
class ConstrainXYFcn : public CharmBaseFcn
{
  public:
    /*
     * Calls parent constructor
     */
    ConstrainXYFcn(const std::vector<double>& data,
                   const std::vector<double>& times,
                   const std::vector<double>& errors,
                   const IntegralOptions_t&   integralOptions);

    /*
     * Creates a model given a vector of {x, y, r, z_im, z_re, width};
     *
     * Returns chi squared value between the model given our params and the measured data, constraining X and Y to the
     * values
     */
    virtual double operator()(const std::vector<double>& parameters) const override;
};

/*
 * Fit ratio of D->K3Pi rates to the equation involving r_D, x, y, etc.; without constraining x and y to their world
 * average values
 */
class ConstrainedFitter : public CharmFitterBase
{
  public:
    /*
     * Tell the fitter about the bin limits and what initial guesses of our parameters and errors to make
     *
     * Parameter ordering is defined by _paramNames
     */
    ConstrainedFitter(const std::vector<double>&   binLimits,
                      const std::array<double, 6>& initialValues,
                      const std::array<double, 6>& initialErrors);

    /*
     * Fit our data to (r^2 + r(yRZ + xImZ)*Gamma*t + (x2 +y2)/4 (Gamma*t)2), constraining x and y to their world
     * average values
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

#endif // CHARM_CONSTRAINED_FITTER_H
