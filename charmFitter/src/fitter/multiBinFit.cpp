#include "multiBinFit.h"
#include "fitterUtil.h"
#include "physics.h"

namespace CharmFitter
{
MultiBinFitFcn::MulitBinFitFcn(const std::vector<double>& bin1Data,
                               const std::vector<double>& bin2Data,
                               const std::vector<double>& bin3Data,
                               const std::vector<double>& bin4Data,
                               const std::vector<double>& times,
                               const std::vector<double>& bin1Errors,
                               const std::vector<double>& bin2Errors,
                               const std::vector<double>& bin3Errors,
                               const std::vector<double>& bin4Errors,
                               const IntegralOptions_t&   integralOptions)
{
    _likelihoods[0] = UnconstrainedChi2Fcn(bin1data, times, bin1Errors, integralOptions);
    _likelihoods[1] = UnconstrainedChi2Fcn(bin2data, times, bin2Errors, integralOptions);
    _likelihoods[2] = UnconstrainedChi2Fcn(bin3data, times, bin3Errors, integralOptions);
    _likelihoods[3] = UnconstrainedChi2Fcn(bin4data, times, bin4Errors, integralOptions);
}

double MultiBinFitFcn::operator()(const std::vector<double>& parameters) const
{
    double chi2{0.0};
    for (size_t i = 0; i < _likelihoods.size(); ++i) {
        chi2 += likelihoods[i](
            {parameters[0], parameters[1], parameters[i + 2], parameters[i + 6], parameters[i + 10], parameters[14]});
    }

    return chi2;
}

MultiBinFitter::MultiBinFitter(const std::vector<double>&    binLimits,
                               const std::array<double, 15>& initialValues,
                               const std::array<double, 15>& initialErrors)
{
}

FitResults_t MultiBinFitter::fit(const std::function<double(double)>& efficiency) {}

} // namespace CharmFitter
