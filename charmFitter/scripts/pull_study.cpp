/*
 * I think the fitter is currently unbiased in the limit of infinite events
 */
#include <iostream>

#include "DecaySimulator.h"
#include "FitterUtils.h"
#include "MinuitPolynomialFitter.h"
#include "PhysicalFitter.h"
#include "RatioCalculator.h"
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
 * Plot a fit
 *
 * This should be two functions
 */
void plotFit(std::array<double, 3>& expectedFitParams,
             const PhysicalFitter&  MyFitter,
             const double           maxTime,
             const size_t           i,
             const bool             drawTrueFit = true)
{
    const util::LegendParams_t legend = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = ""};
    MyFitter.plot->SetTitle("Example Fit with Simulated Data;time/ps;DCS/CF ratio");

    if (drawTrueFit) {
        TF1* trueFit = new TF1("true fit", "[0] +[1]*x+[2]*x*x", 0, maxTime);
        trueFit->SetParameter(0, expectedFitParams[0]);
        trueFit->SetParameter(1, expectedFitParams[1]);
        trueFit->SetParameter(2, expectedFitParams[2]);
        trueFit->SetLineColor(kGray);

        util::saveObjectsToFile<TGraph>(
            std::vector<TObject*>{MyFitter.plot.get(), MyFitter.bestFitFunction.get(), trueFit},
            std::vector<std::string>{"AP", "SAME", "SAME"},
            std::vector<std::string>{"Data", "best fit", "'True' fit"},
            std::string{"fit" + std::to_string(i) + ".pdf"},
            legend);
        delete trueFit;
    } else {
        util::saveObjectsToFile<TGraph>(std::vector<TObject*>{MyFitter.plot.get(), MyFitter.bestFitFunction.get()},
                                        std::vector<std::string>{"AP", "SAME"},
                                        std::vector<std::string>{"Data", "best fit"},
                                        std::string{"fit" + std::to_string(i) + ".pdf"},
                                        legend);
    }
}

/*
 * Generate a vector of {{x vals}, {y vals}} to use in our experiments
 */
std::vector<std::vector<double>> generateXYvals(const std::shared_ptr<std::mt19937>& rndGen,
                                                const size_t                         numExperiments)
{
    std::vector<std::vector<double>> xyCovariance{
        std::vector<double>{WORLD_AVERAGE_X_ERR * WORLD_AVERAGE_X_ERR,
                            X_Y_CORRELATION * WORLD_AVERAGE_X_ERR * WORLD_AVERAGE_Y_ERR},
        std::vector<double>{X_Y_CORRELATION * WORLD_AVERAGE_X_ERR * WORLD_AVERAGE_Y_ERR,
                            WORLD_AVERAGE_Y_ERR * WORLD_AVERAGE_Y_ERR}};

    return util::correlatedGaussianNumbers(
        rndGen, numExperiments, std::vector<double>{WORLD_AVERAGE_X, WORLD_AVERAGE_Y}, xyCovariance);
}

void pull_study(const size_t meanNumCfEvents, const size_t numExperiments)
{
    // Choose parameters to use when simulating
    double width               = 2.4390;
    double maxTime             = 10 / width;
    size_t numBins             = 25;
    double efficiencyTimescale = 1 / width;

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
        double meanNumDcsEvents = Phys::numDCSDecays(meanNumCfEvents, phaseSpaceParams, maxTime, efficiencyTimescale);

        // Generator and PDF for random numbers
        std::uniform_real_distribution<double> uniform;
        auto cfRate  = [&](double x) { return Phys::cfRate(x, phaseSpaceParams, efficiencyTimescale); };
        auto dcsRate = [&](double x) { return Phys::dcsRate(x, phaseSpaceParams, efficiencyTimescale); };

        auto gen = [&](void) {
            double x = uniform(*rndGen);
            double z = 1 - std::exp(-1 * phaseSpaceParams.width * maxTime);
            return (-1 / phaseSpaceParams.width) * std::log(1 - z * x);
        };
        auto genPDF = [&](double x) {
            return std::exp(-phaseSpaceParams.width * x) * phaseSpaceParams.width /
                   (1 - std::exp(-phaseSpaceParams.width * maxTime));
        };
        SimulatedDecays MyDecays = SimulatedDecays(gen, genPDF, cfRate, dcsRate, std::make_pair(0., maxTime), rndGen);

        // Find how many decays to simulate
        std::poisson_distribution<size_t> cfDist(meanNumCfEvents);
        std::poisson_distribution<size_t> dcsDist(meanNumDcsEvents);
        size_t                            numCfEvents  = cfDist(*rndGen);
        size_t                            numDcsEvents = dcsDist(*rndGen);

        // Simulate them
        MyDecays.findCfDecayTimes(numCfEvents);
        MyDecays.findDcsDecayTimes(numDcsEvents);
        MyDecays.plotRates(binLimits);

        // Time binning
        std::vector<size_t> cfCounts  = util::binVector(MyDecays.RSDecayTimes, binLimits);
        std::vector<size_t> dcsCounts = util::binVector(MyDecays.WSDecayTimes, binLimits);

        // Find the ratio
        std::pair<std::vector<double>, std::vector<double>> ratiosAndErrors =
            RatioCalculator::ratioAndError(cfCounts, dcsCounts);

        IntegralOptions_t integralOptions(true, phaseSpaceParams.width, binLimits, 0.5 * efficiencyTimescale);
        FitData_t         MyFitData(binLimits, ratiosAndErrors.first, ratiosAndErrors.second);

        // Fit data
        PhysicalFitter MyFitter(MyFitData, integralOptions, true);
        MyFitter.setFitParams(std::vector<double>{phaseSpaceParams.x,
                                                  phaseSpaceParams.y,
                                                  phaseSpaceParams.r,
                                                  phaseSpaceParams.z_im,
                                                  phaseSpaceParams.z_re,
                                                  phaseSpaceParams.width},
                              std::vector<double>(6, 1));

        // First perform a fit with fixed r, then store the value of chi squared.
        MyFitter.fixParameters(std::vector<std::string>{"width", "z_im", "r"});
        MyFitter.fit();
        double chiSqFixedR = MyFitter.fitParams.fitStatistic;
        MyFitter.freeParameters(std::vector<std::string>{"r"});

        // Now perform a fit with fixed Z and store chisq
        MyFitter.fixParameters(std::vector<std::string>{"z_re"});
        MyFitter.fit();
        double chiSqFixedZ = MyFitter.fitParams.fitStatistic;
        MyFitter.freeParameters(std::vector<std::string>{"z_re"});

        // Finally perform a fit with r and Re(Z) free
        MyFitter.fit();

        // Store the distance of our fit chisq from chisq when r was fixed
        double deltaChiSqR = chiSqFixedR - MyFitter.fitParams.fitStatistic;
        rCoverage[i]       = std::sqrt(std::fabs(deltaChiSqR));
        double deltaChiSqZ = chiSqFixedZ - MyFitter.fitParams.fitStatistic;
        zCoverage[i]       = std::sqrt(std::fabs(deltaChiSqZ));

        auto expectedFitParams = Phys::expectedParams(phaseSpaceParams);
        plotFit(expectedFitParams, MyFitter, maxTime, i);

        // Store parameter and chi squared
        reZPull[i] = (MyFitter.fitParams.fitParams[4] - phaseSpaceParams.z_re) / MyFitter.fitParams.fitParamErrors[4];
        rPull[i]   = (MyFitter.fitParams.fitParams[2] - phaseSpaceParams.r) / MyFitter.fitParams.fitParamErrors[2];
        chiSquaredVals[i] = MyFitter.fitParams.fitStatistic;

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
    pull_study(20000000, 10000);
}
