#include <iostream>
#include <vector>

#include "ReadRoot.h"
#include "efficiencyUtil.h"

#include <TFile.h>
#include <TTree.h>

ReadRoot::ReadRoot(TFile *                         myTFile,
                   const std::string &             treeName,
                   const std::vector<std::string> &branchNames,
                   const std::vector<std::string> &momentumPostfixes)
    : _file(myTFile), _momentumPostfixes(momentumPostfixes)
{
    // Check that our ROOT file exists
    if (!myTFile || myTFile->IsZombie()) {
        throw InvalidRootFile();
    }

    // Check that we have been given 4 branch names
    if (branchNames.size() != 4) {
        throw InvalidBranchesException(branchNames.size());
    }

    // Find what our branch names should be
    std::string kBranchName   = branchNames[0];
    std::string pi1BranchName = branchNames[1];
    std::string pi2BranchName = branchNames[2];
    std::string pi3BranchName = branchNames[3];

    // Allocate memory to our tree (i assume)
    _file->GetObject(treeName.c_str(), _tree);
    if (!_tree) {
        throw InvalidTree(treeName);
    }
    _numEvents = _tree->GetEntries();

    // Initialise vector of D->K3pi events to the right length
    events = std::vector<dDecay_t>(_numEvents);

    // Fill in elements
    // Create buffers for storing our the values in our ROOT file
    // These will be Lorentz vectors for each particle
    std::vector<double> kBuffer(4);
    std::vector<double> pi1Buffer(4);
    std::vector<double> pi2Buffer(4);
    std::vector<double> pi3Buffer(4);
    _setBranchAddresses(kBranchName, kBuffer);
    _setBranchAddresses(pi1BranchName, pi1Buffer);
    _setBranchAddresses(pi2BranchName, pi2Buffer);
    _setBranchAddresses(pi3BranchName, pi3Buffer);

    for (size_t i = 0; i < _numEvents; ++i) {
        // Store our kinematic parameters into the buffers
        _tree->GetEntry(i);

        // Copy from buffers into the right entries of the events vector
        events[i].kParams.px     = kBuffer[0];
        events[i].kParams.py     = kBuffer[1];
        events[i].kParams.pz     = kBuffer[2];
        events[i].kParams.energy = kBuffer[3];

        events[i].pi1Params.px     = pi1Buffer[0];
        events[i].pi1Params.py     = pi1Buffer[1];
        events[i].pi1Params.pz     = pi1Buffer[2];
        events[i].pi1Params.energy = pi1Buffer[3];

        events[i].pi2Params.px     = pi2Buffer[0];
        events[i].pi2Params.py     = pi2Buffer[1];
        events[i].pi2Params.pz     = pi2Buffer[2];
        events[i].pi2Params.energy = pi2Buffer[3];

        events[i].pi3Params.px     = pi3Buffer[0];
        events[i].pi3Params.py     = pi3Buffer[1];
        events[i].pi3Params.pz     = pi3Buffer[2];
        events[i].pi3Params.energy = pi3Buffer[3];
    }
    _tree->ResetBranchAddresses();
}

void ReadRoot::_setBranchAddresses(const std::string &branchPrefix, std::vector<double> &buffer)
{
    // This HAS to be 4 as there are only 4 elements in a momentum 4-vector
    for (size_t i = 0; i < 4; ++i) {
        // Check that our branch exists
        if (!_tree->GetBranch((branchPrefix + _momentumPostfixes[i]).c_str())) {
            throw BranchNotFoundException(branchPrefix + _momentumPostfixes[i]);
        }

        // Make our TBranch point at the right element in the buffer
        _tree->SetBranchAddress((branchPrefix + _momentumPostfixes[i]).c_str(), &buffer[i]);
    }
}
