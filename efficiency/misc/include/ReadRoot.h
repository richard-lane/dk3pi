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
    BranchNotFoundException(const std::string &name) : msg("Failed to find branch " + name) { ; }

    // Should probably also say where we failed to find the branch...
    const char *what() const throw() { return (msg).c_str(); }

    const std::string msg;
};

/*
 * Invalid ROOT file
 */
struct InvalidRootFile : public std::exception {
    InvalidRootFile(void) { ; }

    // The TFile constructor will print a message warning us which file didn't open, so we don't need to also print it
    // here
    const char *what() const throw() { return "Root file failed to open"; }
};

/*
 * Invalid tree
 */
struct InvalidTree : public std::exception {
    InvalidTree(const std::string &treeName) : msg(treeName + " not found") { ; }

    const char *      what() const throw() { return msg.c_str(); }
    const std::string msg;
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
 * Class representing the kinematic data stored in a ROOT file for a series of D->k3pi decays which should really be a
 * function
 *
 * Assumes the data is stored in branches called things like k_Px, k_Py, ... etc.
 * The momentum postfixes (_Px, _Py, _Pz, _E) can be provided to the constructor
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
    explicit ReadRoot(
        TFile *                         myTFile,
        const std::string &             treeName,
        const std::vector<std::string> &branchNames,
        const std::vector<std::string> &momentumPostfixes = std::vector<std::string>{"_Px", "_Py", "_Pz", "_E"});

    /*
     * Vector of events
     */
    std::vector<dDecay_t> events{};

  private:
    /*
     * Take in a branch prefix (might be e.g. _1_K~ for branches called _1_K~_Px, _1_K~_Py, ...) and assign the (px, py,
     * pz, E) parts of the this branch to point to the elements of buffer
     *
     * This is how reading data from a ROOT file works: we give ROOT a location and a branch name and in future it will
     * write data from that branch to that location. In this case, our branches (*_px, *_py, *_pz, *_E) will write to
     * the vector (px, py, pz, E)
     *
     * Maybe this is an inefficient way of doing things but it is ok for now; I/O shouldn't be the slow part of my
     * analysis...
     */
    void _setBranchAddresses(const std::string &branchPrefix, std::vector<double> &buffer);

    TTree *                        _tree;
    TFile *                        _file;
    size_t                         _numEvents;
    const std::vector<std::string> _momentumPostfixes;
};

#endif // READROOT_H