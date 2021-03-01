#ifndef CLEO_INTERFACE_H
#define CLEO_INTERFACE_H

#include "cleo_ll.h"
#include "fitterUtil.h"

#include <map>

namespace CLEO
{

/*
 * CLEO likelihood requires us to provide 64 parameters
 */
constexpr short numParams{64};

/*
 * The CLEO likelihood depends on a number of parameters, whose default values are provided here
 */
constexpr double defaultCleoParams[numParams]{20418.9,    // D0{K-,pi+,pi+,pi-}_D0{K+,pi-,pi-,pi+}_N
                                              17920.6,    // D0{K-,pi+,pi+,pi-}_D0{K+,pi-}_N
                                              4189,       // D0{K-,pi+}_D0{K+,pi-}_N
                                              54259,      // D0{K-,pi+,pi0}_D0{K+,pi-,pi0}_N
                                              63703.4,    // D0{K-,pi+,pi+,pi-}_D0{K+,pi-,pi0}_N
                                              31044.6,    // D0{K-,pi+,pi0}_D0{K+,pi-}_N
                                              0.0553431,  // r_k3pi
                                              0.0440767,  // r_kpipi0
                                              0.627289,   // R_k3pi_0
                                              99.5899,    // d_k3pi_0
                                              0.817555,   // R_kpipi0
                                              202.984,    // d_kpipi0
                                              0.0590658,  // r_kpi
                                              1,          // R_kpi
                                              190.356,    // d_kpi
                                              0.999999,   // R_k3pi_1
                                              132.766,    // d_k3pi_1
                                              0.483193,   // R_k3pi_2
                                              165.398,    // d_k3pi_2
                                              0.260427,   // R_k3pi_3
                                              331.317,    // d_k3pi_3
                                              0.00143751, // D0{pi+,pi-}_BR
                                              0.00681448, // y
                                              0.0790972,  // D0{K-,pi+,pi+,pi-}_BR
                                              0.00406694, // D0{K+,K-}_BR
                                              3249.77,    // D0{K0S0,pi0}_N
                                              0.0119,     // D0{K0S0,pi0}_BR
                                              3467.6,     // D0{K0S0,omega(782)0{pi+,pi-,pi0}}_N
                                              0.0099,     // D0{K0S0,omega(782)0{pi+,pi-,pi0}}_BR
                                              3732.31,    // D0{K0S0,pi0,pi0}_N
                                              0.0091,     // D0{K0S0,pi0,pi0}_BR
                                              4203.02,    // D0{K0S0,phi(1020)0{K+,K-}}_N
                                              0.002,      // D0{K0S0,phi(1020)0{K+,K-}}_BR
                                              3523.35,    // D0{K0S0,eta0}_N
                                              0.00189,    // D0{K0S0,eta0}_BR
                                              2671.25,    // D0{K0S0,eta0{pi+,pi-,pi0}}_N
                                              0.0011,     // D0{K0S0,eta0{pi+,pi-,pi0}}_BR
                                              3082.67,    // D0{K0S0,eta'(958)0{eta0,pi+,pi-}}_N
                                              0.00159,    // D0{K0S0,eta'(958)0{eta0,pi+,pi-}}_BR
                                              3464.94,    // D0{K0L0,pi0}_N
                                              0.01,       // D0{K0L0,pi0}_BR
                                              3469.86,    // D0{K0L0,omega(782)0{pi+,pi-,pi0}}_N
                                              0.0099,     // D0{K0L0,omega(782)0{pi+,pi-,pi0}}_BR
                                              3500.22,    // D0{pi+,pi-,pi0}_N
                                              0.0147,     // D0{pi+,pi-,pi0}_BR
                                              0.973,      // f_pipipi0
                                              2889.28,    // D0{K-,pi+,pi0}_D0{K0S0,pi+,pi-}_N
                                              0.655,      // D0{K0S0,pi+,pi-}::c1
                                              -0.025,     // D0{K0S0,pi+,pi-}::s1
                                              2174.13,    // D0{K0S0,pi+,pi-}_N
                                              0.511,      // D0{K0S0,pi+,pi-}::c2
                                              0.141,      // D0{K0S0,pi+,pi-}::s2
                                              0.024,      // D0{K0S0,pi+,pi-}::c3
                                              1.111,      // D0{K0S0,pi+,pi-}::s3
                                              -0.569,     // D0{K0S0,pi+,pi-}::c4
                                              0.328,      // D0{K0S0,pi+,pi-}::s4
                                              -0.903,     // D0{K0S0,pi+,pi-}::c5
                                              -0.181,     // D0{K0S0,pi+,pi-}::s5
                                              -0.616,     // D0{K0S0,pi+,pi-}::c6
                                              -0.52,      // D0{K0S0,pi+,pi-}::s6
                                              0.1,        // D0{K0S0,pi+,pi-}::c7
                                              -1.129,     // D0{K0S0,pi+,pi-}::s7
                                              0.422,      // D0{K0S0,pi+,pi-}::c8
                                              -0.35};     // D0{K0S0,pi+,pi-}::s8

/*
 * The CLEO paper used four phsp bins
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
 * Find the CLEO likelihood
 *
 * Not multiplied by -2
 */
inline double cleoLikelihood(const Bin phspBin, const FitterUtil::DecayParams_t& params)
{
    // R and d (mag and phase of interference parameter) depend on the bin we're considering...
    // These should really be constexpr, somehow
    std::map<const Bin, const short> rIndices{{Bin::Bin1, 8}, {Bin::Bin2, 15}, {Bin::Bin3, 17}, {Bin::Bin4, 19}};
    std::map<const Bin, const short> dIndices{{Bin::Bin1, 9}, {Bin::Bin2, 16}, {Bin::Bin3, 18}, {Bin::Bin4, 20}};

    // Need to convert Re and Im parts of Z to magnitude and phase
    const std::complex z{params.z_re, params.z_im};
    const double       mag   = std::abs(z);
    const double       phase = 180.0 * std::arg(z) / M_PI;

    // Construct an array of the parameter we want to pass to the CLEO likelihood function
    std::array<double, numParams> cleoParams{};
    std::copy(std::begin(defaultCleoParams), std::end(defaultCleoParams), cleoParams.begin());
    cleoParams[22]                = params.y;
    cleoParams[6]                 = params.r;
    cleoParams[rIndices[phspBin]] = mag;
    cleoParams[dIndices[phspBin]] = phase;

    return cleo_ll(cleoParams.data());
}

} // namespace CLEO

#endif // CLEO_INTERFACE_H
