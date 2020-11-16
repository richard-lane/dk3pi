#ifndef S_WEIGHTING_H
#define S_WEIGHTING_H

#include <RooAbsPdf.h>
#include <RooDstD0BG.h>
#include <RooJohnson.h>
#include <TTree.h>

namespace sWeighting
{

/*
 * Details for a branch describing an observable
 *
 * Stores data as doubles
 */
typedef struct Observable {
    Observable(const std::string& name, const double min, const double max, const std::string& units)
        : name(name), range(std::make_pair(min, max)), units(units)
    {
    }
    const std::string                           name{};
    const std::pair<const double, const double> range{};
    const std::string                           units{};
} Observable_t;

/*
 * Take in an nTuple containing raw (signal + background) D->K3pi data (inFile)
 *
 * Write an nTuple containing (only) the sWeighting branches to outFile
 *
 * Fit signal and background models to it and add sWeight branches
 * Signal and background weights are on branches named numSignalEvents_sw and numBackgroundEvents_sw
 * Likelihoods are on branches names L_numSignalEvents and L_numBackgroundEvents
 *
 * Probably would be a good idea to have a backup of your ROOT file in case this breaks it somehow
 *
 * Pass in the name of the observable parameter in the models as observable. e.g. DELTA_M
 *
 * if a mass fit plot C-string is provided then a plot of the mass fit will be created
 *
 * if a graphViz diagram C-string is provided then a graphviz .dot file showing the model structure will be created
 *
 * multithreaded massfit uses 8 CPUS by default
 *
 */
void findSWeights(const std::string&              inFile,
                  const std::string&              outFile,
                  const std::string&              treeName,
                  RooAbsPdf&                      signalModel,
                  RooAbsPdf&                      backgroundModel,
                  const int                       expectedNumSignal,
                  const int                       expectedNumBackground,
                  const Observable_t&             observable,
                  const std::vector<std::string>& fixedParams,
                  const char*                     massFitPlot     = nullptr,
                  const char*                     graphVizDiagram = nullptr,
                  const int                       numCPU          = 8);

} // namespace sWeighting

#endif // S_WEIGHTING_H
