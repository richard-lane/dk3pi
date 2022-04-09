#include "BesCombinationFitter.h"
#include "fitterUtil.h"
#include "physics.h"

namespace CharmFitter
{
CombinedBESFcn::CombinedBESFcn(const std::vector<double>& data,
                               const std::vector<double>& times,
                               const std::vector<double>& errors,
                               const IntegralOptions_t&   integralOptions,
                               const BES::Bin             binNumber)
    : ConstrainXYFcn(data, times, errors, integralOptions), _binNumber(binNumber)
{
    ;
}

double CombinedBESFcn::operator()(const std::vector<double>& parameters) const
{

    // The BES likelihood may return NaN if Z is outside of the allowed region
    // Assign a special value in this case
    const FitterUtil::DecayParams_t params{
        .x     = parameters[0],
        .y     = parameters[1],
        .r     = parameters[2],
        .z_im  = parameters[3],
        .z_re  = parameters[4],
        .width = parameters[5],
    };
    double besConstraint = BES::besLikelihood(_binNumber, params);
    besConstraint        = std::isnan(besConstraint) ? _nonsenseChi2 : besConstraint;

    // Now find the likelihood based on the D decay times and the world average X and Y
    double chi2 = ConstrainXYFcn::operator()(parameters);

    return chi2 + besConstraint;
}

const double CombinedBESFcn::_nonsenseChi2{10000.0};

BESCombinationFitter::BESCombinationFitter(const std::vector<double>&   binLimits,
                                           const std::array<double, 6>& initialValues,
                                           const std::array<double, 6>& initialErrors,
                                           const BES::Bin               binNumber)
    : ConstrainedFitter(binLimits, initialValues, initialErrors), _binNumber(binNumber)
{
    // Copy the param names, initial values and errors into vectors for passing into _setParams
    const std::vector<double>      vals(initialValues.begin(), initialValues.end());
    const std::vector<double>      errs(initialErrors.begin(), initialErrors.end());
    const std::vector<std::string> names(_paramNames.begin(), _paramNames.end());

    _setParams(names, vals, errs);
}

BESCombinationFitter::BESCombinationFitter(const std::vector<double>&   binLimits,
                                           const std::array<double, 6>& initialValues,
                                           const std::array<double, 6>& initialErrors,
                                           const int                    binNumber)
    : BESCombinationFitter(binLimits, initialValues, initialErrors, static_cast<BES::Bin>(binNumber))
{
}

FitResults_t BESCombinationFitter::fit(const std::function<double(double)>& efficiency)
{
    // This code is exactly the same as for the constrained fitter....
    // TODO refactor to make it nicer
    _fitFcn = std::make_unique<CombinedBESFcn>(
        this->ratios(), this->getBinCentres(), this->errors(), IntegralOptions_t{_binEdges, efficiency}, _binNumber);

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

} // namespace CharmFitter
