#ifndef S_WEIGHTING_H
#define S_WEIGHTING_H

#include <RooAbsPdf.h>
#include <RooDstD0BG.h>
#include <RooJohnson.h>
#include <TTree.h>

namespace sWeighting
{

struct BranchMismatch : public std::exception {
    const char* what() const throw() { return "Number of ROOT branches, allowed ranges and units did not match."; };
};

/*
 * Tell us the name, units, and min/max allowed values of a root branch
 *
 * Only deals with doubles for now
 */
struct RootBranch {
    std::string name{};
    std::string units{};
    double      min{};
    double      max{};
};

/*
 * Take in an nTuple containing raw (signal + background) D->K3pi data (inFile)
 *
 * Fit it to signal and background models, and return a TTree containing the weights and likelihoods
 * Signal and background weights are on branches named numSignalEvents_sw and numBackgroundEvents_sw
 * Likelihoods are on branches names L_numSignalEvents and L_numBackgroundEvents
 *
 * actually it returns a tree containing everything and might break if you give it a large dataset, so. fix that
 *
 * Pass in a vector of (branch name, numerical limits, units) as branchesLimitsAndUnits
 *
 * Pass in the name of the observable parameter in the models as observable. e.g. DELTA_M
 *
 * if a mass fit plot C-string is provided then a plot of the mass fit will be created
 * if a graphViz diagram C-string is provided then a graph showing the model structure will be created
 *
 * multithreaded massfit uses 8 CPUS by default
 *
 */
std::unique_ptr<TTree> findSWeights(const std::string&              inFile,
                                    const std::string&              treeName,
                                    std::vector<RootBranch>         branchesLimitsAndUnits,
                                    RooAbsPdf&                      signalModel,
                                    RooAbsPdf&                      backgroundModel,
                                    const std::string&              observable,
                                    const std::vector<std::string>& fixedParams,
                                    const char*                     massFitPlot     = nullptr,
                                    const char*                     graphVizDiagram = nullptr,
                                    const int                       numCPU          = 8);

} // namespace sWeighting

#endif // S_WEIGHTING_H
