/*
 * Read in the D2K3Pi data from a ROOT file generated by AmpGen
 */
#ifndef READAMPGEN_H
#define READAMPGEN_H

#include "TLorentzVector.h"
#include "TTree.h"

/*
 * Class representing the data stored in a ROOT file for a series of decays
 *
 */
class ReadRoot
{
  public:
    /*
     * All the information we need to extract is in the ROOT file on a given Tree, so this is all the constructor needs.
     *
     * Sets numEvents.
     */
    explicit ReadRoot(TFile *myTFile, std::string treeName);

    size_t numEvents;

  protected:
    TTree *myTree;
};

/*
 * Class for reading D->K3Pi data from AmpGen
 */
class D2K3PiData : public ReadRoot
{
    using ReadRoot::ReadRoot;

  public:
    // Particle data
    std::vector<double>         decayTimes{};
    std::vector<TLorentzVector> kVectors{};
    std::vector<TLorentzVector> pi1Vectors{};
    std::vector<TLorentzVector> pi2Vectors{};
    std::vector<TLorentzVector> pi3Vectors{};

    /*
     * Populate a D2K3PiData class with particle data
     *
     */
    void populate(std::string timesBranchName);

  private:
    /*
     * Write the data on branchName to the index'th position of each TLorentzVector in myVector.
     *
     * e.g. to write x-momenta of a particle described by ROOT branch foo_Px, call writeBranchToLorentzVectors("foo_Px",
     * myVector, 0)
     *
     * The TLorentzVector should be of the form (Px, Py, Pz, E).
     */
    void writeBranchToLorentzVectors(const std::string &          branchName,
                                     std::vector<TLorentzVector> &myVector,
                                     const size_t &               index);
    /*
     * Write the data for a given particle to a vector of TLorentzVectors
     * Relies on particle data being represented by a branch named branchName = "_Px" etc.
     *
     */
    const std::vector<TLorentzVector> particleData(std::string particleName);

    /*
     * Set decay times on a ROOT file branch
     * timesBranchName might be set to e.g. "D_decayTime" for the decay of a D meson
     *
     */
    void setDecayTimes(std::string timesBranchName);
};

#endif // READAMPGEN_H