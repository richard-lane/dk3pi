#include "k3pi_binning.h"

void example_k3pi_binning()
{
    gRandom = new TRandom3(5);
    TGenPhaseSpace phsp;
    TLorentzVector pD(0, 0, 0, 1.86962);
    /// expected ordering in event is K+,pi-,pi-,pi+
    std::vector<double> masses = {0.493677, 0.13957018, 0.13957018, 0.13957018};
    phsp.SetDecay(pD, 4, masses.data());

    // apply scaling and rotation to the DCS amplitude such that we get dcs/cf amplitude ratio 'r' = 0.055
    // and the average relative strong-phase between the two amplitudes ~ 0.
    // 1.04827 0.0601387
    std::complex<double>  dcs_offset = 0.0601387 * exp(std::complex<double>(0, 1) * 1.04827 * M_PI / 180.);
    k3pi_binning::binning bins("binning/dcs.so",   /// dcs model file
                               "binning/cf.so",    /// cf  model file
                               dcs_offset,         /// (constant) offset to DCS amplitude
                               {-39, 0, 43, 180}); /// Bin limits in phase (centred on zero by construction)

    /// calculate global hadronic parameters, and parameters in each of the bins.

    std::complex<double>              z(0, 0);
    double                            n_cf(0);
    double                            n_dcs(0);
    std::vector<std::complex<double>> z_binned(5, std::complex<double>(0, 0));
    std::vector<double>               n_cf_binned(5, 0);
    std::vector<double>               n_dcs_binned(5, 0);
    for (int i = 0; i < 5 * pow(10, 6); ++i) {
        auto event_4v = k3pi_binning::makeUnweighted(phsp);
        auto event    = k3pi_binning::eventFromVectors(event_4v);
        auto eval_cf  = bins.cf(event.data(), 1);
        auto eval_dcs = dcs_offset * bins.dcs(event.data(), 1);
        z += eval_cf * std::conj(eval_dcs);
        n_cf += std::norm(eval_cf);
        n_dcs += std::norm(eval_dcs);

        auto bin = bins.bin(event_4v, 1);

        z_binned[bin] += eval_cf * std::conj(eval_dcs);
        n_cf_binned[bin] += std::norm(eval_cf);
        n_dcs_binned[bin] += std::norm(eval_dcs);
    }
    std::cout << "==== Global: =====================" << std::endl;
    std::cout << "R = " << std::abs(z) / sqrt(n_cf * n_dcs) << std::endl; /// should be ~ 0.47
    std::cout << "d = " << std::arg(z) * 180 / M_PI << std::endl; /// should be ~ zero by construction (< 1 degree)
    std::cout << "r = " << sqrt(n_dcs / n_cf) << std::endl;       /// should be ~ 0.055

    for (int i = 0; i < 5; ++i) {
        std::cout << "==== Bin " << i + 1 << ": ======================" << std::endl;
        std::cout << "R[" << i + 1 << "] = " << std::abs(z_binned[i]) / sqrt(n_cf_binned[i] * n_dcs_binned[i])
                  << std::endl;
        std::cout << "d[" << i + 1 << "] = " << std::arg(z_binned[i]) * 180 / M_PI << std::endl;
        std::cout << "K[" << i + 1 << "] = " << n_dcs_binned[i] / n_dcs << std::endl;
        std::cout << "K'[" << i + 1 << "] = " << n_cf_binned[i] / n_cf << std::endl;
    }
    std::cout << "==================================" << std::endl;
}
