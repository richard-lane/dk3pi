#include <iostream>
#include <vector>

#include "DecaySimulator.h"
#include "Fitter.h"
#include "PullStudyHelpers.h"
#include "RatioCalculator.h"
#include "util.h"

#include "TGraph.h"

void splitVectorOfPairs(std::vector<std::pair<double, double>>& pairs,
                        std::vector<double>&                    first,
                        std::vector<double>&                    second)
{
    for (auto it = std::make_move_iterator(pairs.begin()), end = std::make_move_iterator(pairs.end()); it != end;
         ++it) {
        first.push_back(std::move(it->first));
        second.push_back(std::move(it->second));
    }
}

void plotScan(Fitter& MinuitChiSqFitter, const size_t numPoints, const size_t paramIndex)
{
    std::string graphTitle{""};
    switch (paramIndex) {
    case 0: graphTitle = "Parameter a ChiSq scan"; break;
    case 1: graphTitle = "Parameter b ChiSq scan"; break;
    case 2: graphTitle = "Parameter c ChiSq scan"; break;
    default: throw;
    }

    MinuitChiSqFitter.chiSqParameterScan(paramIndex, numPoints);
    std::vector<double> values{};
    std::vector<double> chiSq{};
    splitVectorOfPairs(MinuitChiSqFitter.parameterScan, values, chiSq);
    TGraph* Graph = new TGraph(numPoints, values.data(), chiSq.data());
    Graph->SetTitle((graphTitle + ";value;chiSq").c_str());
    util::saveObjectToFile(Graph, (graphTitle + ".pdf").c_str(), "AP");
}

/*
 * Create a dataset using accept-reject and our RatioCalculator then fit it to polynomials using
 * Minuit ChiSq
 *
 * Perform a scan of each parameter, store them as std::vector<std::pair<double, double>> and plot
 */
void test_param_scan(void)
{
    // Create an accept-reject dataset
    DecayParams_t phaseSpaceParams = {
        .x     = 0.0037,
        .y     = 0.0066,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = 2439.0,
    };
    double                    maxTime      = 0.002;
    std::pair<double, double> allowedTimes = std::make_pair(0, maxTime);
    std::pair<double, double> allowedRates = std::make_pair(0, 1.3);

    size_t numCfEvents  = 1000000;
    size_t numDcsEvents = PullStudyHelpers::numDCSDecays(numCfEvents, phaseSpaceParams, maxTime);

    SimulatedDecays MyDecays = SimulatedDecays(allowedTimes, allowedRates, phaseSpaceParams);
    MyDecays.findDcsDecayTimes(numDcsEvents);
    MyDecays.findCfDecayTimes(numCfEvents);

    // Define some time bins
    std::vector<double> dcsTimes{MyDecays.WSDecayTimes};
    std::sort(dcsTimes.begin(), dcsTimes.end());
    std::vector<double> timeBinLimits = util::findBinLimits(dcsTimes, 100, 0, 1.05 * maxTime);

    // Divide using RatioCalculator
    RatioCalculator MyRatios = RatioCalculator(MyDecays.RSDecayTimes, MyDecays.WSDecayTimes, timeBinLimits);
    MyRatios.calculateRatios();

    // Create a fitter
    FitData_t MyFitData         = FitData(MyRatios.binCentres, MyRatios.binWidths, MyRatios.ratio, MyRatios.error);
    Fitter    MinuitChiSqFitter = Fitter(MyFitData);

    // Perform fit, outputu minimum statistic
    std::vector<double> initialParameterGuess{0.02, 1.0, 100.0};
    std::vector<double> initialErrorsGuess{0.01, 1.0, 100.0};
    MinuitChiSqFitter.fitUsingMinuit(initialParameterGuess, initialErrorsGuess, ChiSquared);
    std::cout << "Min chisq: " << *(MinuitChiSqFitter.statistic) << std::endl;
    std::cout << "Params: " << MinuitChiSqFitter.fitParams.fitParams[0] << " "
              << MinuitChiSqFitter.fitParams.fitParams[1] << " " << MinuitChiSqFitter.fitParams.fitParams[2]
              << std::endl;

    // Perform a chi squared scan on each parameter
    size_t numPoints = 100;
    plotScan(MinuitChiSqFitter, numPoints, 0);
    plotScan(MinuitChiSqFitter, numPoints, 1);
    plotScan(MinuitChiSqFitter, numPoints, 2);
}

int main()
{
    test_param_scan();
    return 0;
}
