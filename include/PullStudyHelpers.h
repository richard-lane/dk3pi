#ifndef PULL_STUDY_HELPERS_H
#define PULL_STUDY_HELPERS_H

#include <string>
#include <vector>

#include "DecaySimulator.h"

/*
 * Plot the distribution of a vector as a histogram of 200 bins
 *
 * Plot is centred on expectedMean and spans 10*expectedSigma
 * Sometimes it doesn't though because root is odd
 */
void plot_parameter_distribution(std::string         title,
                                 std::vector<double> parameter,
                                 size_t              nExperiments,
                                 double              expectedMean  = 0,
                                 double              expectedSigma = 1);

/*
 * Find our expected a, b and c
 */
std::vector<double> expectedParams(const DecayParams_t &phaseSpaceParams);

/*
 * Find the number of DCS decays we need to simulate, given the number of CF decays and our phase space parameters.
 *
 * The ratio of DCS to CF decays is calculated from the ratio of the integrals of their decay rates.
 *
 */
size_t numDCSDecays(const size_t numCFDecays, const DecayParams_t &phaseSpaceParams, double maxTime);

#endif // PULL_STUDY_HELPERS_H
