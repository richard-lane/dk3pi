#ifndef CUTS_H
#define CUTS_H
#include <string>

/*
 * Copy the tree from inFile to outFile, subject to some cuts
 *
 * Look at the implementation for the cuts cus they're subject to change
 */
void promptCuts(const std::string& inFile,
                const std::string& treeName,
                const std::string& outFile,
                const std::string& piSoftPTBranchName,
                const std::string& piSoftGhostProbBranchName,
                const std::string& dPtBranchName,
                const std::string& dMassBranchName);

/*
 * Copy the tree from inFile to outFile, subject to some cuts
 *
 * Look at the implementation for the cuts cus they're subject to change
 */
void semileptonicCuts(const std::string& inFile,
                      const std::string& treeName,
                      const std::string& outFile,
                      const std::string& dPtBranchName,
                      const std::string& dMassBranchName);

#endif // CUTS_H
