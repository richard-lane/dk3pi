#ifndef CHARM_POLY_FITTER_H
#define CHARM_POLY_FITTER_H

#include "CharmFitterBase.h"

namespace CharmFitter
{

/*
 * Minuit2-compatible class for fitting to a polynomial using chi squared
 *
 * Fits by integrating over the bins
 */
class PolynomialChiSqFcn : public CharmBaseFcn
{
  public:
    /*
     * Calls parent constructor
     *
     * Must also provide a D decay width to this fitter, as this isn't a fit parameter and we need it to estimate the WS
     * and RS rates.
     */
    PolynomialChiSqFcn(const std::vector<double>& data,
                       const std::vector<double>& times,
                       const std::vector<double>& errors,
                       const double               width,
                       const IntegralOptions_t&   integralOptions);

    /*
     * Creates a model given a vector of three parameters.
     *
     * Returns the chi squared value between the model given our parameters and the measured data.
     */
    virtual double operator()(const std::vector<double>& parameters) const override;

  private:
    double _width{};
};

/*
 * Fit ratio of D->K3Pi rates to a + bt + ct^2
 */
class CharmPolynomialFitter : public CharmFitterBase
{
  public:
    /*
     * Tell the fitter about the bin limits and what initial guesses of our parameters and errors to make
     *
     * Must also provide a D decay width, since this is needed by the fitter
     */
    CharmPolynomialFitter(const std::vector<double>&   binLimits,
                          const std::array<double, 3>& initialValues,
                          const std::array<double, 3>& initialErrors,
                          const double                 width);
    /*
     * Fit our data to a second-order polynomial a + bt + ct^2 using Minuit2 and the chi-squared method.
     *
     */
    FitResults_t fit(const std::function<double(double)>& efficiency) override;

  private:
    double _width{};
};

} // namespace CharmFitter
#endif // CHARM_POLY_FITTER_H
