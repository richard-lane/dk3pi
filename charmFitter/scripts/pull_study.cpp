/*
 * I think the fitter is currently unbiased in the limit of infinite events
 */
#include <iostream>

#include "ConstrainedFitter.h"
#include "DecaySimulator.h"
#include "FitterUtils.h"
#include "PolynomialFitter.h"
#include "RatioCalculator.h"
#include "UnconstrainedFitter.h"
#include "physics.h"
#include "util.h"

#include <TH1D.h>
#include <TMath.h>

#include <boost/progress.hpp>

void plot_parameter_distribution(std::string         title,
                                 std::vector<double> parameter,
                                 size_t              nExperiments,
                                 double              expectedMean,
                                 double              expectedSigma)
{
    // Define axis limits
    double xMin = expectedMean - 5 * expectedSigma;
    double xMax = expectedMean - 5 * expectedSigma;

    TH1D* MyGraph = new TH1D(title.c_str(), title.c_str(), 200, xMin, xMax);

    MyGraph->FillN(nExperiments, parameter.data(), 0);
    MyGraph->SetTitle((title + ";Normalised Pull;Count").c_str());

    util::saveObjectToFile(MyGraph, (title + ".pdf").c_str());

    std::cout << title + " mean:\t\t" + MyGraph->GetMean() << std::endl;
    std::cout << title + " std dev:\t" + MyGraph->GetStdDev() << std::endl;
    delete MyGraph;
}

/*
 * Generate a vector of {{x vals}, {y vals}} to use in our experiments
 */
std::vector<std::vector<double>> generateXYvals(const std::shared_ptr<std::mt19937>& rndGen,
                                                const size_t                         numExperiments)
{
    std::vector<std::vector<double>> xyCovariance{
        std::vector<double>{CharmFitter::WORLD_AVERAGE_X_ERR * CharmFitter::WORLD_AVERAGE_X_ERR,
                            CharmFitter::X_Y_CORRELATION * CharmFitter::WORLD_AVERAGE_X_ERR *
                                CharmFitter::WORLD_AVERAGE_Y_ERR},
        std::vector<double>{CharmFitter::X_Y_CORRELATION * CharmFitter::WORLD_AVERAGE_X_ERR *
                                CharmFitter::WORLD_AVERAGE_Y_ERR,
                            CharmFitter::WORLD_AVERAGE_Y_ERR * CharmFitter::WORLD_AVERAGE_Y_ERR}};

    return util::correlatedGaussianNumbers(
        rndGen,
        numExperiments,
        std::vector<double>{CharmFitter::WORLD_AVERAGE_X, CharmFitter::WORLD_AVERAGE_Y},
        xyCovariance);
}

void pull_study(const size_t meanNumCfEvents, const size_t numExperiments)
{
    // Choose parameters to use when simulating
    double width   = 2.4390;
    double maxTime = 10 / width;
    size_t numBins = 25; // not actually the number of bins wtf

    // Create RNGs for numbers of decays
    std::random_device            rd;
    std::shared_ptr<std::mt19937> rndGen = std::make_shared<std::mt19937>(rd());

    // Find exponentially-spaced time bin limits to use
    // Make some bins at the start wider (because of the efficiency)
    std::vector<double> binLimits = FitterUtil::exponentialBinLimits(maxTime, width, numBins);
    binLimits.erase(binLimits.begin() + 1, binLimits.begin() + 3);
    binLimits.erase(binLimits.begin() + 4);

    // Generate X and Y values to be used in the simulations
    std::vector<std::vector<double>> xyVals = generateXYvals(rndGen, numExperiments);

    // Initialise vectors of fit parameter pulls and chi squared
    std::vector<double> rPull(numExperiments, -1);
    std::vector<double> reZPull(numExperiments, -1);
    std::vector<double> chiSquaredVals(numExperiments, -1);

    // Initialise vector of delta chi^2 values to test the r distribution
    std::vector<double> rCoverage(numExperiments, -1);
    std::vector<double> zCoverage(numExperiments, -1);

    // No efficiency function
    auto efficiency = [](const double) { return 1; };

    boost::progress_display showProgress(numExperiments);

    for (size_t i = 0; i < numExperiments; ++i) {
        double                    thisX            = xyVals[0][i];
        double                    thisY            = xyVals[1][i];
        FitterUtil::DecayParams_t phaseSpaceParams = {
            .x     = thisX,
            .y     = thisY,
            .r     = 0.055,
            .z_im  = -0.2956,
            .z_re  = 0.7609,
            .width = width,
        };
        double meanNumDcsEvents = Phys::numDCSDecays(meanNumCfEvents, phaseSpaceParams, maxTime);

        // Generator and PDF for random numbers
        SimulatedDecays MyDecays = SimulatedDecays({0, maxTime}, phaseSpaceParams, *rndGen);

        // Find how many decays to simulate
        std::poisson_distribution<size_t> cfDist(meanNumCfEvents);
        std::poisson_distribution<size_t> dcsDist(meanNumDcsEvents);
        size_t                            numCfEvents  = cfDist(*rndGen);
        size_t                            numDcsEvents = dcsDist(*rndGen);

        // Create a fitter here
        // Have to do it inside the loop as the phaseSpaceParams were generated randomly
        CharmFitter::ConstrainedFitter MyFitter(binLimits,
                                                {phaseSpaceParams.x,
                                                 phaseSpaceParams.y,
                                                 phaseSpaceParams.r,
                                                 phaseSpaceParams.z_im,
                                                 phaseSpaceParams.z_re,
                                                 phaseSpaceParams.width},
                                                {1, 1, 1, 1, 1, 1});

        // Tell the fitter about the data
        MyFitter.addRSPoints(MyDecays.rsDecayTimes(numCfEvents), std::vector<double>(numCfEvents, 1.0));
        MyFitter.addWSPoints(MyDecays.wsDecayTimes(numDcsEvents), std::vector<double>(numDcsEvents, 1.0));

        // First perform a fit with fixed r, then store the value of chi squared.
        MyFitter.fixParameters(std::array{"width", "z_im", "r"});
        const double chiSqFixedR = MyFitter.fit(efficiency).fitStatistic;
        MyFitter.freeParameter("r");

        // Now perform a fit with fixed Z and store chisq
        MyFitter.fixParameter("z_re");
        const double chiSqFixedZ = MyFitter.fit(efficiency).fitStatistic;
        MyFitter.freeParameter("z_re");

        // Finally perform a fit with r and Re(Z) free
        // Store the distance of our fit chisq from chisq when r was fixed
        const auto&  result      = MyFitter.fit(efficiency);
        const double deltaChiSqR = chiSqFixedR - result.fitStatistic;
        rCoverage[i]             = std::sqrt(std::fabs(deltaChiSqR));
        const double deltaChiSqZ = chiSqFixedZ - result.fitStatistic;
        zCoverage[i]             = std::sqrt(std::fabs(deltaChiSqZ));

        // Store parameter and chi squared
        reZPull[i]        = (result.fitParams[4] - phaseSpaceParams.z_re) / result.fitParamErrors[4];
        rPull[i]          = (result.fitParams[2] - phaseSpaceParams.r) / result.fitParamErrors[2];
        chiSquaredVals[i] = result.fitStatistic;

        ++showProgress;
    }

    // Output mean and std dev of pulls
    std::vector<std::pair<double, double>> stats{};

    stats.push_back(util::meanAndStdDev(rPull));
    stats.push_back(util::meanAndStdDev(reZPull));

    std::cout << "Constrained XY Fit:" << std::endl;
    for (auto pair = stats.begin(); pair != stats.end(); ++pair) {
        std::cout << "Pull:\t" << pair->first << "+-" << pair->second << std::endl;
    }
    std::cout << std::endl;

    plot_parameter_distribution("r", rPull, numExperiments, 0, 1);
    plot_parameter_distribution("Re(Z)", reZPull, numExperiments, 0, 1);

    util::saveHistogram(chiSquaredVals, "MinuitChiSq.png", "", 100);

    // Plot a cumulative histogram of delta chi squared values
    size_t              nBins{50};
    double              max{4};
    std::vector<double> bins(nBins + 1);
    for (size_t i = 0; i <= nBins; ++i) {
        bins[i] = i * (max / nBins);
    }

    TH1* hr = new TH1D("rCumulative", "r coverage", nBins, bins.data());
    hr->FillN(numExperiments, rCoverage.data(), nullptr);
    TH1* hrc = hr->GetCumulative();
    hrc->SetStats(false);

    TH1* hz = new TH1D("zCumulative", "Re(Z) coverage", nBins, bins.data());
    hz->FillN(numExperiments, zCoverage.data(), nullptr);
    TH1* hzc = hz->GetCumulative();
    hzc->SetStats(false);

    TF1* expectedErf = new TF1("erf", "(TMath::Freq(x) - TMath::Freq(0))*2*[0]", 0, max);
    expectedErf->SetParameter(0, numExperiments);

    const util::LegendParams_t legend = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = ""};
    util::saveObjectsToFile<TGraph>(std::vector<TObject*>{hrc, expectedErf},
                                    std::vector<std::string>{"", "CSAME"},
                                    std::vector<std::string>{"delta Chisq", "expected"},
                                    "rCoverage.pdf",
                                    legend);

    util::saveObjectsToFile<TGraph>(std::vector<TObject*>{hzc, expectedErf},
                                    std::vector<std::string>{"", "CSAME"},
                                    std::vector<std::string>{"delta Chisq", "expected"},
                                    "zCoverage.pdf",
                                    legend);
}

int main()
{
    pull_study(5000000, 100);
}
