#include <TF1.h>

#include "UnconstrainedFitter.h"
#include "fitterUtil.h"
#include "physics.h"

namespace CharmFitter
{

UnconstrainedChiSqFcn::UnconstrainedChiSqFcn(const std::vector<double>& data,
                                             const std::vector<double>& times,
                                             const std::vector<double>& errors,
                                             const IntegralOptions_t&   integralOptions)
    : CharmBaseFcn(data, times, errors, integralOptions)
{
    ;
}

double UnconstrainedChiSqFcn::operator()(const std::vector<double>& parameters) const
{
    // Generate lambdas for our RS and WS rates, given this set of parameters
    auto rsRate = [&parameters](const double x) { return Phys::cfRate(x, parameters[5]); };
    auto wsRate = [&parameters](const double x) { return Phys::dcsRate(x, parameters); };

    double chi2 = 0.0;
    for (size_t i = 0; i < theMeasurements.size(); ++i) {
        double model =
            util::gaussLegendreQuad(wsRate, _integralOptions.binLimits[i], _integralOptions.binLimits[i + 1]) /
            util::gaussLegendreQuad(rsRate, _integralOptions.binLimits[i], _integralOptions.binLimits[i + 1]);
        chi2 += std::pow((model - theMeasurements[i]) / theMVariances[i], 2);
    }

    return chi2;
}

UnconstrainedFitter::UnconstrainedFitter(const std::vector<double>&   binLimits,
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

FitResults_t UnconstrainedFitter::fit(const std::function<double(double)>& integralOptions)
{
    _fitFcn = std::make_unique<UnconstrainedChiSqFcn>(
        this->ratios(), this->getBinCentres(), this->errors(), IntegralOptions_t{_binEdges, integralOptions});

    // Check that we have fixed at least three of x, y, Re(Z), Im(Z) and decay width.
    // Copy the parameters we have into a vector + remove r
    std::vector<ROOT::Minuit2::MinuitParameter> paramsExceptR = _parameters->Trafo().Parameters();
    paramsExceptR.erase(paramsExceptR.begin() + 2);
    // NB: r being fixed does not help us, as we lose 1 piece of information (the constant term)

    // Count how many of the relevant parameters are fixed
    unsigned short numFixedParams{0};
    for (auto it = paramsExceptR.begin(); it != paramsExceptR.end(); ++it) {
        if (it->IsFixed()) {
            ++numFixedParams;
        }
    }

    // Raise if not enough are fixed
    if (numFixedParams < 3) {
        std::cerr << "Must fix at least three of {x, y, width, Im(Z), Re(Z)} for fit to be well defined ("
                  << numFixedParams << " fixed)." << std::endl;
        throw D2K3PiException();
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

const std::array<std::string, 6> UnconstrainedFitter::_paramNames{"x", "y", "r", "z_im", "z_re", "width"};

} // namespace CharmFitter
