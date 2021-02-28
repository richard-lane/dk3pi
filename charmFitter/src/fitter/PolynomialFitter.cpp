#include <TF1.h>

#include "PolynomialFitter.h"
#include "physics.h"

namespace CharmFitter
{

PolynomialChiSqFcn::PolynomialChiSqFcn(const std::vector<double>& data,
                                       const std::vector<double>& times,
                                       const std::vector<double>& errors,
                                       const double               width,
                                       const IntegralOptions_t&   integralOptions)
    : CharmBaseFcn(data, times, errors, integralOptions), _width(width)
{
    ;
}

double PolynomialChiSqFcn::operator()(const std::vector<double>& parameters) const
{
    // Generate lambdas for our RS and WS rates, given this set of parameters
    auto rsRate = [this](const double x) { return Phys::cfRate(x, _width); };
    auto wsRate = [&parameters, this](const double x) { return Phys::dcsRatePoly(x, parameters, _width); };

    double chi2 = 0.0;
    for (size_t i = 0; i < theMeasurements.size(); ++i) {
        double model =
            util::gaussLegendreQuad(wsRate, _integralOptions.binLimits[i], _integralOptions.binLimits[i + 1]) /
            util::gaussLegendreQuad(rsRate, _integralOptions.binLimits[i], _integralOptions.binLimits[i + 1]);
        chi2 += std::pow((model - theMeasurements[i]) / theMVariances[i], 2);
    }
    return chi2;
}

CharmPolynomialFitter::CharmPolynomialFitter(const std::vector<double>&   binLimits,
                                             const std::array<double, 3>& initialValues,
                                             const std::array<double, 3>& initialErrors,
                                             const double                 width)
    : CharmFitterBase(binLimits), _width(width)
{
    // Copy the initial values and errors into a vector for passing into _setParams
    const std::vector<double> vals(initialValues.begin(), initialValues.end());
    const std::vector<double> errs(initialErrors.begin(), initialErrors.end());

    _setParams({"a", "b", "c"}, vals, errs);
}

FitResults_t CharmPolynomialFitter::fit(const std::function<double(double)>& efficiency)
{
    _fitFcn = std::make_unique<PolynomialChiSqFcn>(
        this->ratios(), this->getBinCentres(), this->errors(), _width, IntegralOptions_t{_binEdges, efficiency});
    FitResults_t result = CharmFitterBase::_performFit();

    result.bestFitFunction = TF1("Best fit function", "[0] +[1]*x+[2]*x*x");
    result.bestFitFunction.SetParameters(result.fitParams.data());

    return result;
}

} // namespace CharmFitter
