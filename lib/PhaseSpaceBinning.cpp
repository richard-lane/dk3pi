#ifndef PHASE_SPACE_BINNING_CPP
#define PHASE_SPACE_BINNING_CPP

#include <boost/filesystem.hpp>
#include <iostream>
#include <vector>

#include "TLorentzVector.h"

#include "boost/progress.hpp"

#include "D2K3PiError.h"
#include "PhaseSpaceBinning.h"
#include "k3pi_binning.h"

PhaseSpaceBinning::PhaseSpaceBinning(const std::vector<double> &        decayTimes,
                                     const std::vector<TLorentzVector> &kVectors,
                                     const std::vector<TLorentzVector> &pi1Vectors,
                                     const std::vector<TLorentzVector> &pi2Vectors,
                                     const std::vector<TLorentzVector> &pi3Vectors)
{
    _decayTimes = decayTimes;
    _kVectors   = kVectors;
    _pi1Vectors = pi1Vectors;
    _pi2Vectors = pi2Vectors;
    _pi3Vectors = pi3Vectors;

    // TODO add error checking for different length vectors
    _numEvents = decayTimes.size();
}

void PhaseSpaceBinning::performBinning()
{
    // Check that _dcsFile and _cfFile exist
    if (!boost::filesystem::exists(_dcsFile) || !boost::filesystem::exists(_cfFile)) {
        std::cerr << _dcsFile << " or " << _cfFile << " not found." << std::endl;
        throw D2K3PiException();
    }

    // Progress bar
    boost::progress_display progressBar(_numEvents);

    // Define our bins
    k3pi_binning::binning bins(_dcsFile, _cfFile, _dcs_offset, {BIN_LIMITS});

    // Iterate over all events, sorting them into bins
    for (size_t i = 0; i < _numEvents; ++i) {
        // Create a vector of TLorentzVectors for this event (K+, pi-, pi-, pi+)
        std::vector<TLorentzVector> eventVector{_kVectors[i], _pi1Vectors[i], _pi2Vectors[i], _pi3Vectors[i]};
        auto                        event = k3pi_binning::eventFromVectors(eventVector);

        // Find which bin the event belongs in
        // The 1 tags the sign of the K meson in the D->K3pi decay
        int bin = bins.bin(eventVector, 1);

        // Log the time of the event in this bin
        // This might be slow because it's dynamically resizing the vector, but should be ok for our purposes.
        // TODO A marginally cleverer implementation might set binnedTimes[bin] to the number of points we have, and
        // then just reshape the whole thing once we've finished binning.
        _binnedTimes[bin].push_back(_decayTimes[i]);

        ++progressBar;
    }

    // Output how many points are in each bin
    for (size_t bin = 0; bin < NUM_BINS; ++bin) {
        std::cout << "points in bin " << bin << ": " << _binnedTimes[bin].size() << std::endl;
    }
}

#endif // PHASE_SPACE_BINNING_CPP
