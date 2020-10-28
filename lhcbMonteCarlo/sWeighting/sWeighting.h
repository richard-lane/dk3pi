#ifndef S_WEIGHTING_H
#define S_WEIGHTING_H

#include <RooAbsPdf.h>
#include <RooDstD0BG.h>
#include <RooJohnson.h>
#include <TTree.h>

namespace sWeighting
{

/*
 * Johnson function to use as signal model
 *
 * TODO maybe ive named these parameters wrong
 *
 *     dMassRange: variable encoding the low/high cutoffs for D mass
 *     meanDMass:  the central value of the D mass peak
 *     dMassWidth: width param for D mass peak; corresponds to parameter lambda in the Johnson function
 *     gamma: Johnson function gamma parameter; controls the shape of the distribution, distorts it left/right
 *     delta: Johnson function delta parameter; controls the strength of the Gaussian contribution
 *     massThreshhold: mass below which the PDF is set to 0
 */
RooJohnson johnsonSignalModel(RooAbsReal& dMassRange,
                              RooAbsReal& meanDMass,
                              RooAbsReal& dMassWidth,
                              RooAbsReal& gamma,
                              RooAbsReal& delta,
                              double      massThreshhold);

/*
 * Prompt D background model
 *
 * This is some ridiculous thing that looks like
 *     (1-exp(-(deltaM - deltaM0)/c))*(deltaM/deltaM0)^a  + (deltaM/deltaM0- 1)*b
 *
 * deltaM is the Dstar-D0 mass difference; deltaM0 is the threshhold mass difference- i.e. the low mass difference
 * cutoff
 *
 * a, b, c are shape parameters
 */
RooDstD0BG promptBackgroundModel(RooAbsReal& deltaM, RooAbsReal& deltaM0, RooAbsReal& a, RooAbsReal& b, RooAbsReal& c);

/*
 * If no branch containing the D mass difference exists, we may have to create one
 *
 * This does that
 *
 * Copies the old root file to copy_<filename>, in case of things going bad
 */
void addDeltaMBranch(const std::string& path,
                     const std::string& treeName,
                     const std::string& dMassBranchName,
                     const std::string& dStarMassBranchName,
                     const std::string& deltaMBranchName);

/*
 * Take in an nTuple containing raw (signal + background) D->K3pi data (inFile)
 *
 * Fit it to signal and background models, and write a new nTuple containing the sWeighted tree
 *
 * Should really take in a collection of branch names to preserve, but this is hard-coded for now
 *
 * I think the signal and background models both need to contain a variable named "DELTA_M", and the ROOT file needs a
 * branch named DELTA_M
 * but i don't actually know
 */
void createWeightedRootFile(const std::string&              inFile,
                            const std::string&              treeName,
                            RooAbsPdf&                      signalModel,
                            RooAbsPdf&                      backgroundModel,
                            const std::vector<std::string>& fixedParams,
                            const std::string&              outTreePath);

} // namespace sWeighting

#endif // S_WEIGHTING_H
