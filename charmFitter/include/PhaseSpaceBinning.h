/*
 * Bin D->K3Pi decays in phase space using kinematic data.
 */
#ifndef PHASE_SPACE_BIN_HPP
#define PHASE_SPACE_BIN_HPP

// ---- Magic Numbers
// DCS and CF relative amplitude and phase
#define DCS_MAGNITUDE 0.0445
#define DCS_PHASE -3.04

// Bin limits in phase, centred on zero by construction
#define NUM_BINS 5
#define BIN_LIMITS -39, 0, 43, 180

#include <complex>
#include <string>
#include <vector>

#include "TLorentzVector.h"

/*
 * Bin D->K3Pi decays in phase space using decay time and kinematic particle data
 * Uses Tim Evans' k3pi_binning.cpp code.
 */
class PhaseSpaceBinning
{
  public:
    /*
     * Set the kinematic data that will be used for binning
     */
    PhaseSpaceBinning(const std::vector<double>&         decayTimes,
                      const std::vector<TLorentzVector>& kVectors,
                      const std::vector<TLorentzVector>& pi1Vectors,
                      const std::vector<TLorentzVector>& pi2Vectors,
                      const std::vector<TLorentzVector>& pi3Vectors);
    /*
     * Bin events into predefined phase space bins.
     */
    void performBinning();

    /*
     * Vector of vectors holding binned times.
     */
    std::vector<std::vector<double>> _binnedTimes{NUM_BINS};

  private:
    // Particle data
    std::vector<double>         _decayTimes{};
    std::vector<TLorentzVector> _kVectors{};
    std::vector<TLorentzVector> _pi1Vectors{};
    std::vector<TLorentzVector> _pi2Vectors{};
    std::vector<TLorentzVector> _pi3Vectors{};

    // Parameters needed for binning
    const std::string          _dcsFile{"../AmpGenTools/dcs.so"};
    const std::string          _cfFile{"../AmpGenTools/cf.so"};
    const std::complex<double> _dcs_offset = DCS_MAGNITUDE * exp(std::complex<double>(0, 1) * DCS_PHASE * M_PI / 180.);
    size_t                     _numEvents{0};
};

#endif // PHASE_SPACE_BIN_HPP