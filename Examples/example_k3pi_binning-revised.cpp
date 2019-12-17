/*
 * example_k3pi_binning-revised.cpp
 * An attempt to make sense of how the binning scheme works
 * The same functionality as example_k3pi_binning.cpp, but maybe a bit more readable
 */

#include "k3pi_binning.h"

// ---- Magic numbers
// Decay Product data
#define NUMBER_PRODUCTS 4
#define K_MASS 0.493677
#define PION_MASS 0.13957018

// Decaying particle (Momentum, Energy) vector in GeV
#define PARTICLE_MOMENTUM 0, 0, 0, 1.86962

// DCS and CF relative amplitude and phase
#define DCS_MAGNITUDE 0.0601387
#define DCS_PHASE 1.04827

/// Bin limits in phase, centred on zero by construction
#define NUM_BINS 5
#define BIN_LIMITS -39, 0, 43, 180

// Number of events
#define NUM_EVENTS 5 * pow(10, 6)

void example_k3pi_binning_revised()
{
    TGenPhaseSpace phsp;
    TLorentzVector pD(PARTICLE_MOMENTUM);

    /// expected ordering in event is K+,pi-,pi-,pi+
    std::vector<double> masses = {K_MASS, PION_MASS, PION_MASS, PION_MASS};
    phsp.SetDecay(pD, NUMBER_PRODUCTS, masses.data());

    // apply scaling and rotation to the DCS amplitude such that we get dcs/cf amplitude ratio 'r' = 0.055
    // and the average relative strong-phase between the two amplitudes ~ 0.
    std::complex<double>  dcs_offset = DCS_MAGNITUDE * exp(std::complex<double>(0, 1) * DCS_PHASE * M_PI / 180.);
    k3pi_binning::binning bins("binning/dcs.so", /// dcs model file
                               "binning/cf.so",  /// cf  model file
                               dcs_offset,       /// (constant) offset to DCS amplitude
                               {BIN_LIMITS});

    /// calculate global hadronic parameters, and parameters in each of the bins.
    std::complex<double>              z(0, 0);
    double                            n_cf(0);
    double                            n_dcs(0);
    std::vector<std::complex<double>> z_binned(NUM_BINS, std::complex<double>(0, 0));
    std::vector<double>               n_cf_binned(NUM_BINS, 0);
    std::vector<double>               n_dcs_binned(NUM_BINS, 0);

    for (int i = 0; i < NUM_EVENTS; ++i) {
        // Generate an event and convert it to the type needed for the binning by using eventFromVectors()
        auto event_4v = k3pi_binning::makeUnweighted(phsp);
        auto event    = k3pi_binning::eventFromVectors(event_4v);

        //@@@ work out the CF and DCS amplitudes of the ith event
        // Work out the CF and DCS amplitudes of this event
        auto eval_cf  = bins.cf(event.data(), 1);
        auto eval_dcs = dcs_offset * bins.dcs(event.data(), 1);

        // Update the global hadronic parameters
        z += eval_cf * std::conj(eval_dcs);
        n_cf += std::norm(eval_cf);
        n_dcs += std::norm(eval_dcs);

        // Find which bin the event belongs in and update its hadronic parameters
        auto bin = bins.bin(event_4v, 1);
        z_binned[bin] += eval_cf * std::conj(eval_dcs);
        n_cf_binned[bin] += std::norm(eval_cf);
        n_dcs_binned[bin] += std::norm(eval_dcs);
    }

    std::cout << "==== Global: =====================" << std::endl;
    std::cout << "R = " << std::abs(z) / sqrt(n_cf * n_dcs) << std::endl; /// should be ~ 0.47
    std::cout << "d = " << std::arg(z) * 180 / M_PI << std::endl; /// should be ~ zero by construction (< 1 degree)
    std::cout << "r = " << sqrt(n_dcs / n_cf) << std::endl;       /// should be ~ 0.055

    for (int i = 0; i < NUM_BINS; ++i) {
        std::cout << "==== Bin " << i + 1 << ": ======================" << std::endl;
        std::cout << "R[" << i + 1 << "] = " << std::abs(z_binned[i]) / sqrt(n_cf_binned[i] * n_dcs_binned[i])
                  << std::endl;
        std::cout << "d[" << i + 1 << "] = " << std::arg(z_binned[i]) * 180 / M_PI << std::endl;
        std::cout << "K[" << i + 1 << "] = " << n_dcs_binned[i] / n_dcs << std::endl;
        std::cout << "K'[" << i + 1 << "] = " << n_cf_binned[i] / n_cf << std::endl;
    }
    std::cout << "==================================" << std::endl;
}
