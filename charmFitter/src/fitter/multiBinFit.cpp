#include "multiBinFit.h"
#include "fitterUtil.h"
#include "physics.h"

namespace CharmFitter
{
MultiBinFitFcn::MultiBinFitFcn(const std::vector<double>& bin1Data,
                               const std::vector<double>& bin2Data,
                               const std::vector<double>& bin3Data,
                               const std::vector<double>& bin4Data,
                               const std::vector<double>& times,
                               const std::vector<double>& bin1Errors,
                               const std::vector<double>& bin2Errors,
                               const std::vector<double>& bin3Errors,
                               const std::vector<double>& bin4Errors,
                               const IntegralOptions_t&   integralOptions)
    : theMeasurements1(bin1Data), theMeasurements2(bin2Data), theMeasurements3(bin3Data), theMeasurements4(bin4Data),
      thePositions(times), theMVariances1(bin1Errors), theMVariances2(bin2Errors), theMVariances3(bin3Errors),
      theMVariances4(bin4Errors), theErrorDef(1.), _integralOptions(integralOptions)
{
    _likelihoods.emplace_back(UnconstrainedChiSqFcn(bin1Data, times, bin1Errors, _integralOptions));
    _likelihoods.emplace_back(UnconstrainedChiSqFcn(bin2Data, times, bin2Errors, _integralOptions));
    _likelihoods.emplace_back(UnconstrainedChiSqFcn(bin3Data, times, bin3Errors, _integralOptions));
    _likelihoods.emplace_back(UnconstrainedChiSqFcn(bin4Data, times, bin4Errors, _integralOptions));
}

double MultiBinFitFcn::operator()(const std::vector<double>& parameters) const
{
    double chi2{0.0};
    for (size_t i = 0; i < _likelihoods.size(); ++i) {
        chi2 += _likelihoods[i](
            {parameters[0], parameters[1], parameters[i + 2], parameters[i + 6], parameters[i + 10], parameters[14]});
    }

    return chi2;
}

MultiBinFitter::MultiBinFitter([[maybe_unused]] const std::vector<double>&    binLimits,
                               [[maybe_unused]] const std::array<double, 15>& initialValues,
                               [[maybe_unused]] const std::array<double, 15>& initialErrors)
{
}

FitResults_t MultiBinFitter::fit([[maybe_unused]] const std::function<double(double)>& efficiency)
{
    return FitResults_t{};
}

} // namespace CharmFitter
