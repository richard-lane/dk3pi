#include <TF1.h>

#include "ConstrainedFitter.h"
#include "fitterUtil.h"
#include "physics.h"

namespace CharmFitter
{
ConstrainXYFcn::ConstrainXYFcn(const std::vector<double>& data,
                               const std::vector<double>& times,
                               const std::vector<double>& errors,
                               const IntegralOptions_t&   integralOptions)
    : UnconstrainedChiSqFcn(data, times, errors, integralOptions)
{
    ;
}

double ConstrainXYFcn::operator()(const std::vector<double>& parameters) const
{
    double chi2 = UnconstrainedChiSqFcn::operator()(parameters);

    // Introduce constraint by modifying chi squared
    double dx         = parameters[0] - WORLD_AVERAGE_X;
    double dy         = parameters[1] - WORLD_AVERAGE_Y;
    double constraint = (1 / (1 - X_Y_CORRELATION * X_Y_CORRELATION)) *
                        (std::pow(dx / WORLD_AVERAGE_X_ERR, 2) + std::pow(dy / WORLD_AVERAGE_Y_ERR, 2) -
                         2 * X_Y_CORRELATION * dx * dy / (WORLD_AVERAGE_Y_ERR * WORLD_AVERAGE_X_ERR));

    return chi2 + constraint;
}

ConstrainedFitter::ConstrainedFitter(const std::vector<double>&   binLimits,
                                     const std::array<double, 6>& initialValues,
                                     const std::array<double, 6>& initialErrors)
    : CharmFitterBase(binLimits)
{
    // Copy the param names, initial values and errors into vectors for passing into _setParams
    const std::vector<double>      vals(initialValues.begin(), initialValues.end());
    const std::vector<double>      errs(initialErrors.begin(), initialErrors.end());
    const std::vector<std::string> names(_paramNames.begin(), _paramNames.end());

    _setParams(names, vals, errs);
}

FitResults_t ConstrainedFitter::fit(const std::function<double(double)>& efficiency)
{
    _fitFcn = std::make_unique<ConstrainXYFcn>(
        this->ratios(), this->getBinCentres(), this->errors(), IntegralOptions_t{_binEdges, efficiency});

    // Check we have fixed at least one of {Re(Z), Im(Z) or width}
    std::vector<ROOT::Minuit2::MinuitParameter> reZ_imZ_width{_parameters->Trafo().Parameters()[3],
                                                              _parameters->Trafo().Parameters()[4],
                                                              _parameters->Trafo().Parameters()[5]};
    if (std::all_of(reZ_imZ_width.begin(), reZ_imZ_width.end(), [](ROOT::Minuit2::MinuitParameter& param) {
            return !param.IsFixed();
        })) {
        std::cerr << "Must fix one component of Z or width for fit to be well defined" << std::endl;
    }

    // Use base class implementation to actually perform the fit
    FitResults_t result = CharmFitterBase::_performFit();

    // Create a best-fit function
    // Its range should cover the range of times that our time bins do
    const double minTime   = _binEdges.front();
    const double maxTime   = _binEdges.back();
    result.bestFitFunction = TF1("Best fit function",
                                 "[2]*[2] + [2]*([1]*[4] + [0]*[3])*[5]*x + 0.25 * ([0]*[0] + [1]*[1])*[5]*[5]*x*x",
                                 minTime,
                                 maxTime);
    result.bestFitFunction.SetParameters(result.fitParams.data());

    return result;
}

const std::array<std::string, 6> ConstrainedFitter::_paramNames{"x", "y", "r", "z_im", "z_re", "width"};

} // namespace CharmFitter
