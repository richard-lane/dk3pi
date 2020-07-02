#ifndef READROOT_H
#define READROOT_H

#include <string>

#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>

#include "efficiencyUtil.h"

/*
 * Couldn't find a branch of the specified name
 */
struct BranchNotFoundException : public std::exception {
    BranchNotFoundException(const std::string &name) : name(name) { ; }

    // Should probably also say where we failed to find the branch...
    const char *what() const throw() { return ("Failed to find branch " + name).c_str(); }

    const std::string name;
};

/*
 * Wrong number of branch names passed in to constructor
 */
struct InvalidBranchesException : public std::exception {
    InvalidBranchesException(const size_t numBranchesProvided) : numBranchesProvided(numBranchesProvided) { ; }

    const char *what() const throw()
    {
        return ("Must provide 4 branches to read K3pi parameters from; received " + (char)numBranchesProvided);
    }
    const size_t numBranchesProvided;
};

/*
 * Class representing the data stored in a ROOT file for a series of D->k3pi decays
 *
 */
class ReadRoot
{
  public:
    /*
     * Pass in a ROOT TFile as well as the name of the tree where the D->k3pi data is stored.
     *
     * Need to pass in a vector of ROOT branch names in the form {kBranchName, pi1, pi2, pi3} for a D -> K pi1 pi2 pi3
     * decay
     * "_Px", "_Py", "_Pz" and "_E" will be appeded to these branch names when looking for branches... i.e. these should
     * be such that the branches <kBranchName>_Px, <kBranchName>_Py etc. exist in the root file
     *
     */
    explicit ReadRoot(TFile *myTFile, const std::string &treeName, const std::vector<std::string> &branchNames);

    /*
     * Vector of events
     */
    std::vector<dDecay_t> events{};

  private:
    /*
     * Take in a branch prefix (e.g. _1_K~ for a K-) and assign the (px, py, pz, E) parts of the this branch to point to
     * the elements of buffer
     */
    void _setBranchAddresses(const std::string &branchPrefix, std::vector<double> &buffer);

    TTree *_tree;
    TFile *_file;
    size_t _numEvents;
};

#endif // READROOT_H