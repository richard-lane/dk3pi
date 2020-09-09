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

PhaseSpaceBinning::PhaseSpaceBinning(const PhspBinningData&     data,
                                     const std::string&         dcsPath,
                                     const std::string          cfPath,
                                     const std::vector<double>  binLimits,
                                     const std::complex<double> dcsOffset)
    : _data(data), _binLimits(binLimits), _dcsModel(dcsPath), _cfModel(cfPath), _dcsOffset(dcsOffset)
{
    ;
}

std::vector<PhspBinningData> PhaseSpaceBinning::performBinning()
{
    // Check that _dcsFile and _cfFile exist
    if (!boost::filesystem::exists(_dcsModel) || !boost::filesystem::exists(_cfModel)) {
        std::cerr << _dcsModel << " or " << _cfModel << " not found." << std::endl;
        throw D2K3PiException();
    }

    // Progress bar
    size_t                  numEvents{_data.decayTimes.size()};
    boost::progress_display progressBar(numEvents);

    // Define our bins
    k3pi_binning::binning bins(_dcsModel, _cfModel, _dcsOffset, _binLimits);

    // Iterate over all events, sorting them into bins
    std::vector<PhspBinningData> binnedData(_numBins);
    for (size_t i = 0; i < numEvents; ++i) {
        // Create a vector of TLorentzVectors for this event (K+, pi-, pi-, pi+)
        std::vector<TLorentzVector> eventVector{
            _data.kVectors[i], _data.pi1Vectors[i], _data.pi2Vectors[i], _data.pi3Vectors[i]};
        auto event = k3pi_binning::eventFromVectors(eventVector);

        // Find which bin the event belongs in
        // The 1 tags the sign of the K meson in the D->K3pi decay
        int bin = bins.bin(eventVector, 1);

        // Log the time of the event in this bin
        // This might be slow because it's dynamically resizing the vector, but should be ok for our purposes.
        // TODO A marginally cleverer implementation might set binnedTimes[bin] to the number of points we have, and
        // then just reshape the whole thing once we've finished binning.
        binnedData[bin].decayTimes.push_back(_data.decayTimes[i]);
        binnedData[bin].kVectors.push_back(_data.kVectors[i]);
        binnedData[bin].pi1Vectors.push_back(_data.pi1Vectors[i]);
        binnedData[bin].pi2Vectors.push_back(_data.pi2Vectors[i]);
        binnedData[bin].pi3Vectors.push_back(_data.pi3Vectors[i]);

        ++progressBar;
    }

    // Output how many points are in each bin
    for (size_t bin = 0; bin < _numBins; ++bin) {
        std::cout << "points in bin " << bin << ": " << binnedData[bin].decayTimes.size() << std::endl;
    }

    return binnedData;
}

#endif // PHASE_SPACE_BINNING_CPP
