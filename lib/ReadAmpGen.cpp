#ifndef READAMPGEN_CPP
#define READAMPGEN_CPP

#include <vector>

#include "../include/ReadAmpGen.h"

#include "TFile.h"

ReadRoot::ReadRoot(TFile *myTFile, std::string treeName)
{
    myTFile->GetObject(treeName.c_str(), myTree);
    numEvents = myTree->GetEntries();
}

void D2K3PiData::populate(std::string timesBranchName)
{
    kVectors   = particleData("_1_K~");
    pi1Vectors = particleData("_2_pi#");
    pi2Vectors = particleData("_3_pi#");
    pi3Vectors = particleData("_4_pi~");
    setDecayTimes(timesBranchName);
}

void D2K3PiData::writeBranchToLorentzVectors(const std::string &          branchName,
                                             std::vector<TLorentzVector> &myVector,
                                             const size_t &               index)
{
    double myData{0.0};

    myTree->SetBranchAddress(branchName.c_str(), &myData);
    for (size_t i = 0; i < numEvents; ++i) {
        myTree->GetEntry(i);
        myVector[i][index] = myData;
    }

    // Reset all branch addresses to avoid a bug where repeatedly calling this function would set an array to the wrong
    // values
    myTree->ResetBranchAddresses();
}

const std::vector<TLorentzVector> D2K3PiData::particleData(std::string particleName)
{
    std::vector<TLorentzVector> myVector = std::vector<TLorentzVector>(numEvents);

    writeBranchToLorentzVectors(particleName + "_Px", myVector, 0);
    writeBranchToLorentzVectors(particleName + "_Py", myVector, 1);
    writeBranchToLorentzVectors(particleName + "_Pz", myVector, 2);
    writeBranchToLorentzVectors(particleName + "_E", myVector, 3);

    return myVector;
}

void D2K3PiData::setDecayTimes(std::string timesBranchName)
{

    // Init our decay times to -1; it should then be obvious if something has gone wrong.
    double myData{0.0};
    decayTimes = std::vector<double>(numEvents, -1);

    myTree->SetBranchAddress(timesBranchName.c_str(), &myData);
    for (size_t i = 0; i < numEvents; ++i) {
        myTree->GetEntry(i);
        decayTimes[i] = myData;
    }

    // Reset all branch addresses to avoid a bug where repeatedly calling this function would set an array to the wrong
    // values
    myTree->ResetBranchAddresses();
}

#endif // READAMPGEN_CPP