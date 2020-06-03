#ifndef PULL_STUDY_HELPERS_H
#define PULL_STUDY_HELPERS_H

#include <string>
#include <vector>

#include "DecaySimulator.h"

namespace PullStudyHelpers
{

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
 * Plot a histogram from a vector
 */
void plotHist(const std::vector<double>& vector, const size_t numBins, const std::string& name);

/*
 * Find the number of DCS decays we need to simulate, given the number of CF decays and our phase space parameters.
 *
 * The ratio of DCS to CF decays is calculated from the ratio of the integrals of their decay rates.
 * CF rate is exponential; DCS is exp * (a + bt + ct^2)
 *
 * Returns a double so e.g. can be used as the mean of a distribution. Cast to an integer type before using as a
 * count!
 *
 */
double numDCSDecays(const size_t         numCFDecays,
                    const DecayParams_t& phaseSpaceParams,
                    double               maxTime,
                    double               efficiencyTimescale);

} // namespace PullStudyHelpers

#endif // PULL_STUDY_HELPERS_H
