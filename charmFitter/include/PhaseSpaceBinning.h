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
#define BIN_LIMITS      \
    {                   \
        -39, 0, 43, 180 \
    }

#include <cassert>
#include <complex>
#include <string>
#include <vector>

#include "TLorentzVector.h"

/*
 * A collection of D->K3Pi events
 *
 * Since we're likely to access all the parameters for a given decay at once,
 * maybe this would be better implemented as a vector of individual events rather than multiple vectors
 *
 * But who knows how TLorentzVector is implemented so it probably doesn't matter
 */
struct PhspBinningData {
    PhspBinningData(void){};
    PhspBinningData(const std::vector<double>&         decayTimes,
                    const std::vector<TLorentzVector>& kVectors,
                    const std::vector<TLorentzVector>& pi1Vectors,
                    const std::vector<TLorentzVector>& pi2Vectors,
                    const std::vector<TLorentzVector>& pi3Vectors)
        : decayTimes(decayTimes), kVectors(kVectors), pi1Vectors(pi1Vectors), pi2Vectors(pi2Vectors),
          pi3Vectors(pi3Vectors)
    {
        size_t numEvents{decayTimes.size()};
        assert(kVectors.size() == numEvents);
        assert(pi1Vectors.size() == numEvents);
        assert(pi2Vectors.size() == numEvents);
        assert(pi3Vectors.size() == numEvents);
    };

    std::vector<double>         decayTimes{};
    std::vector<TLorentzVector> kVectors{};
    std::vector<TLorentzVector> pi1Vectors{};
    std::vector<TLorentzVector> pi2Vectors{};
    std::vector<TLorentzVector> pi3Vectors{};
};

/*
 * Bin D->K3Pi decays in phase space using decay time and kinematic particle data
 * Uses Tim Evans' k3pi_binning.cpp code.
 *
 * Will need to store all our events in memory, which hopefully won't be a problem...
 */
class PhaseSpaceBinning
{
  public:
    /*
     * Set the kinematic and bin data, and tell us where the DCS and CF model files are
     *
     */
    PhaseSpaceBinning(const PhspBinningData&     data,
                      const std::string&         dcsPath,
                      const std::string          cfPath,
                      const std::vector<double>  binLimits = {BIN_LIMITS},
                      const std::complex<double> dcsOffset = std::polar<double>(DCS_MAGNITUDE, DCS_PHASE));

    /*
     * Bin events into predefined phase space bins.
     */
    std::vector<PhspBinningData> performBinning();

  private:
    /*
     * Particle data
     */
    PhspBinningData _data;

    // Parameters needed for binning
    const std::vector<double>  _binLimits;
    size_t                     _numBins = _binLimits.size() + 1;
    const std::string          _dcsModel;
    const std::string          _cfModel;
    const std::complex<double> _dcsOffset;
};

#endif // PHASE_SPACE_BIN_HPP