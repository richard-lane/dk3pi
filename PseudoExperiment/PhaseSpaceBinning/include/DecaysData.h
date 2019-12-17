#ifndef BIN_GENERATED_DECAYS_HPP
#define BIN_GENERATED_DECAYS_HPP

#include <complex>
#include <string>

#include "TLorentzVector.h"
#include "TTree.h"

// ---- Magic Numbers
// DCS and CF relative amplitude and phase
#define DCS_MAGNITUDE 0.0445
#define DCS_PHASE -3.04

// Bin limits in phase, centred on zero by construction
#define NUM_BINS 5
#define BIN_LIMITS -39, 0, 43, 180

/*
 * Class representing the data stored in a ROOT file for a series of decays
 *
 */
class DecaysData
{
  public:
    explicit DecaysData(TFile *myTFile, std::string treeName);
    const std::vector<TLorentzVector> particleData(std::string particleName);

    size_t numEvents;

  protected:
    TTree *myTree;

  private:
    void getNumEvents();

    // Helper for writing vectors
    void writeBranchToLorentzVectors(const std::string &          branchName,
                                     std::vector<TLorentzVector> &myVector,
                                     const size_t &               index);
};

/*
 * Class for D or Dbar decays to K3pi
 */
class D2K3PiData : public DecaysData
{
    using DecaysData::DecaysData;

  public:
    // Particle data
    std::vector<double>         decayTimes{};
    std::vector<TLorentzVector> kVectors{};
    std::vector<TLorentzVector> pi1Vectors{};
    std::vector<TLorentzVector> pi2Vectors{};
    std::vector<TLorentzVector> pi3Vectors{};

    // Phase space binning
    void performBinning(std::string timesBranchName);

    std::vector<std::vector<double>> binnedTimes{NUM_BINS};

    // Time binning
    void plotBinnedTimes(size_t bin);

    std::vector<double>                           timeBinLimits;
    std::vector<std::vector<std::vector<double>>> timeBinnedTimes{NUM_BINS};
    std::vector<std::vector<size_t>>              numPointsPerTimeBin{NUM_BINS};
    std::vector<double>                           binCentres{};
    std::vector<double>                           binErrors{};

  private:
    const std::string          dcsFile{"binning/dcs.so"};
    const std::string          cfFile{"binning/cf.so"};
    const std::complex<double> dcs_offset = DCS_MAGNITUDE * exp(std::complex<double>(0, 1) * DCS_PHASE * M_PI / 180.);

    // Particle data methods
    void populate(std::string timesBranchName);
    void setDecayTimes(std::string timesBranchName);

    // Phase-space binning
    void binTimes(void);

    // Time binning
    size_t numTimeBins{0};

    void                             setTimeBins();
    void                             sortBinnedTimes();
    void                             splitTimes();
    std::vector<std::vector<double>> splitVectorWithLimits(std::vector<double> &myVector);
    void                             setNumPointsPerTimeBin();
    void                             setTimeBinCentres();
};

#endif // BIN_GENERATED_DECAYS_HPP