#ifndef BES_COMBINATION_FITTER_H
#define BES_COMBINATION_FITTER_H

#include "ConstrainedFitter.h"
#include "bes_interface.h"

namespace CharmFitter
{

/*
 * Minuit2-compatible class for fitting charm decay lifetimes, constraining
 * mixing parameters x and y to their world averages and combining with the CLEO charm interference constraint
 *
 * Fits by integrating over the bins
 */
class CombinedBESFcn : public ConstrainXYFcn
{
  public:
    /*
     * Calls parent constructor
     *
     * Must provide the BES phsp bin being considered, as the BES likelihood is evaluated in bins
     */
    CombinedBESFcn(const std::vector<double>& data,
                   const std::vector<double>& times,
                   const std::vector<double>& errors,
                   const IntegralOptions_t&   integralOptions,
                   const BES::Bin             binNumber);

    /*
     * Creates a model given a vector of {x, y, r, z_im, z_re, width};
     *
     * Returns chi squared value between the model given our params and the measured data, constraining X and Y to the
     * world values
     *
     * The chi squared is also constrained based on the BES likelihood
     */
    virtual double operator()(const std::vector<double>& parameters) const override;

  private:
    const BES::Bin _binNumber{};

    /*
     * The chi2 value returned in the event that our BES likelihood evaluates to NaN
     */
    const static double _nonsenseChi2;
};

/*
 * Fit ratio of D->K3Pi rates to the equation involving r_D, x, y, etc.; constraining x and y to their world
 * average values
 */
class BESCombinationFitter : public ConstrainedFitter
{
  public:
    /*
     * Tell the fitter about the bin limits and what initial guesses of our parameters and errors to make
     *
     * Parameter ordering is defined by _paramNames
     *
     * Must provide the BES phsp bin being considered
     */
    BESCombinationFitter(const std::vector<double>&   binLimits,
                         const std::array<double, 6>& initialValues,
                         const std::array<double, 6>& initialErrors,
                         const BES::Bin               binNumber);

    /*
     * Same as above but construct with an int; for the python binding
     */
    BESCombinationFitter(const std::vector<double>&   binLimits,
                         const std::array<double, 6>& initialValues,
                         const std::array<double, 6>& initialErrors,
                         const int                    binNumber);

    /*
     * Fit our data to (r^2 + r(yRZ + xImZ)*Gamma*t + (x2 +y2)/4 (Gamma*t)2), constraining x and y to their world
     * average values
     *
     * Also takes into account the BES likelihood - the BSE likelihood may evaluate to NaN if e.g. some of the
     * parameters are outside of their physically allowed ranges. If this is the case, the BSE contribution to the
     * likelihood is replaced with CombinedBESFcn::_nonsenseChi2
     */
    FitResults_t fit(const std::function<double(double)>& efficiency) override;

  private:
    const BES::Bin _binNumber{};
};

} // namespace CharmFitter

#endif // BES_COMBINATION_FITTER_H
