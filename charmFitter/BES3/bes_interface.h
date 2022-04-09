/*
 * Mostly copy and pasted from the CLEO likelihood
 */

#ifndef BES_INTERFACE_H
#define BES_INTERFACE_H

#include "BESIII_chi2.h"
#include "fitterUtil.h"

#include <map>

namespace BES
{

/*
 * BES 3 likelihood requires us to provide 78 parameters
 */
constexpr short numParams{64};

/*
 * The BES-3 likelihood depends on a number of parameters, whose default values are provided here
 */
constexpr double defaultBESParams[numParams]{
    5.81976e-01,  // Rk3pi0
    7.80911e-01,  // Rk3pi1
    8.49155e-01,  // Rk3pi2
    4.53797e-01,  // Rk3pi3
    1.31418e+02,  // deltak3pi0
    1.49528e+02,  // deltak3pi1
    1.76478e+02,  // deltak3pi2
    2.73767e+02,  // deltak3pi3
    4.41481e-02,  // rkpipi0
    7.95940e-01,  // Rkpipi0
    2.00422e+02,  // deltakpipi0
    9.12111e+03,  // Norm0kspipipi
    7.58751e+03,  // Norm1kspipipi
    7.86921e+03,  // Norm2kspipipi
    1.03944e+04,  // Norm3kspipipi
    1.97604e+02,  // deltakpi
    3.93451e-02,  // Bkpi_CF
    1.47674e-01,  // Bkpipi0_CF
    8.32820e-02,  // Bk3pi_CF
    2.14488e-03,  // Bkpipi0_DCS/Bkpipi0_CF
    3.22223e-03,  // Bk3pi_DCS/Bk3pi_CF
    3.44355e-03,  // rkpi^2
    5.43595e-02,  // rk3pi0
    5.81052e-02,  // rk3pi1
    5.75272e-02,  // rk3pi2
    5.09254e-02,  // rk3pi3
    3.65404e-03,  // x_mixing
    6.87437e-03,  // y_mixing
    6.99649e-01,  // ci_1
    6.45961e-01,  // ci_2
    -2.58903e-03, // ci_3
    -6.07998e-01, // ci_4
    -9.54417e-01, // ci_5
    -5.89476e-01, // ci_6
    5.91171e-02,  // ci_7
    4.11390e-01,  // ci_8
    8.54049e-02,  // si_1
    3.10419e-01,  // si_2
    1.00831e+00,  // si_3
    6.48371e-01,  // si_4
    -3.37068e-02, // si_5
    -5.68738e-01, // si_6
    -8.32666e-01, // si_7
    -4.31009e-01, // si_8
    1.73402e-01,  // ki_1
    8.76005e-02,  // ki_2
    6.91972e-02,  // ki_3
    2.55000e-02,  // ki_4
    8.50011e-02,  // ki_5
    5.91996e-02,  // ki_6
    1.26900e-01,  // ki_7
    1.33800e-01,  // ki_8
    7.93998e-02,  // ki_-1
    1.74000e-02,  // ki_-2
    2.02000e-02,  // ki_-3
    1.62001e-02,  // ki_-4
    5.12002e-02,  // ki_-5
    1.42999e-02,  // ki_-6
    1.32000e-02,  // ki_-7
    2.74999e-02,  // ki_-8
    9.72667e-01,  // Fpipipi0
    6.12508e+04,  // Norm_kspipivskpipi0
    1.75216e-01,  // kip_1
    8.81138e-02,  // kip_2
    6.97188e-02,  // kip_3
    2.54630e-02,  // kip_4
    8.51398e-02,  // kip_5
    5.86190e-02,  // kip_6
    1.26872e-01,  // kip_7
    1.33835e-01,  // kip_8
    8.01789e-02,  // kip_-1
    1.74489e-02,  // kip_-2
    2.01895e-02,  // kip_-3
    1.67866e-02,  // kip_-4
    5.19655e-02,  // kip_-5
    1.32375e-02,  // kip_-6
    1.32823e-02,  // kip_-7
    2.70838e-02   // kip_-8
};

/*
 * The BES paper uses four phsp bins
 *
 * Events are binned based on the phase of the interference parameter
 *
 * Phase should be provided in degrees
 */
enum class Bin { Bin1, Bin2, Bin3, Bin4 };
inline Bin binNumber(const double phase)
{
    // Normalise phase to [-180, 180)
    double normPhase = fmod(phase + 180.0, 360.0);
    if (normPhase < 0)
        normPhase += 360.0;
    normPhase -= 180.0;

    // Find which bin the phase belongs to
    if (normPhase < 39.0)
        return Bin::Bin1;
    if (normPhase < 0.0)
        return Bin::Bin2;
    if (normPhase < 43.0)
        return Bin::Bin3;
    return Bin::Bin4;
}

/*
 * Find the BES likelihood
 *
 * Not multiplied by -2
 */
inline double besLikelihood(const Bin phspBin, const FitterUtil::DecayParams_t& params)
{
    // R and d (mag and phase of interference parameter) depend on the bin we're considering...
    // These should really be constexpr, somehow
    std::map<const Bin, const short> rIndices{{Bin::Bin1, 0}, {Bin::Bin2, 1}, {Bin::Bin3, 2}, {Bin::Bin4, 3}};
    std::map<const Bin, const short> dIndices{{Bin::Bin1, 4}, {Bin::Bin2, 5}, {Bin::Bin3, 6}, {Bin::Bin4, 7}};

    // I guess they measured r_D in each bin?
    std::map<const Bin, const short> rdIndices{{Bin::Bin1, 22}, {Bin::Bin2, 23}, {Bin::Bin3, 24}, {Bin::Bin4, 25}};

    // Need to convert Re and Im parts of Z to magnitude and phase
    const std::complex z{params.z_re, params.z_im};
    const double       mag   = std::abs(z);
    const double       phase = 180.0 * std::arg(z) / M_PI;

    // Construct an array of the parameter we want to pass to the BES likelihood function
    std::array<double, numParams> besParams{};
    std::copy(std::begin(defaultBESParams), std::end(defaultBESParams), besParams.begin());
    besParams[26] = params.x;
    besParams[27] = params.y;

    besParams[rIndices[phspBin]] = mag;
    besParams[dIndices[phspBin]] = phase;

    return BESIII_chi2(besParams.data());
}

} // namespace BES

#endif // BES_INTERFACE_H
