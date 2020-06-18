#ifndef DECAYSIMULATOR_CPP
#define DECAYSIMULATOR_CPP

#include <cmath>
#include <iostream>
#include <utility>

#include <TGraphErrors.h>
#include <TH1D.h>

#include "D2K3PiError.h"
#include "DecaySimulator.h"
#include "physics.h"
#include "util.h"

#include <boost/math/tools/minima.hpp>

SimulatedDecays::SimulatedDecays(const std::function<double(void)> &  generateTime,
                                 const std::function<double(double)> &generatingPDF,
                                 const std::function<double(double)> &cfRate,
                                 const std::function<double(double)> &dcsRate,
                                 const std::pair<double, double> &    timeDomain,
                                 const std::shared_ptr<std::mt19937> &generator)
    : _minTime(timeDomain.first), _maxTime(timeDomain.second)
{
    // Set our generator to the provided random device
    _gen = generator;

    // Set generating functions
    _getRandomTime = generateTime;
    _generatingPDF = generatingPDF;

    // Set rate equations
    _cfRate  = cfRate;
    _dcsRate = dcsRate;

    // Check that our generating pdf is normalised correctly
    if (std::abs(util::gaussLegendreQuad(_generatingPDF, _minTime, _maxTime) - 1.0) > DBL_EPSILON) {
        std::cerr << "generating pdf not normalised; integral between " << _minTime << " and " << _maxTime
                  << " evaluates to " << util::gaussLegendreQuad(_generatingPDF, _minTime, _maxTime) << std::endl;
        throw D2K3PiException();
    }

    // Find the maximum rates; we need these to perform the  accept-reject
    _setMaxRatios();
}

void SimulatedDecays::test(const size_t numPoints, const std::vector<double> &binLimits)
{
    size_t              numBins = binLimits.size() - 1;
    std::vector<double> binCentres(numBins, -1);
    std::vector<double> binWidths(numBins, -1);
    for (size_t i = 0; i < numBins; ++i) {
        binCentres[i] = (binLimits[i] + binLimits[i + 1]) / 2;
        binWidths[i]  = 0.5 * (binLimits[i + 1] - binLimits[i]);
    }

    std::vector<double> expectedNormalisedBinPop(numBins, -1);
    for (size_t i = 0; i < numBins; ++i) {
        expectedNormalisedBinPop[i] = util::gaussLegendreQuad(_generatingPDF, binLimits[i], binLimits[i + 1]);
    }

    // Generate numPoints points and bin them
    std::vector<double> generatedTimes(numPoints, -1);
    for (size_t i = 0; i < numPoints; ++i) {
        generatedTimes[i] = _getRandomTime();
    }
    std::vector<size_t> generatedBinPopulations = util::binVector(generatedTimes, binLimits);

    // Cast to a double cus im tired and dont want to think about it too hard
    std::vector<double> normalisedBinPops(numBins, -1);
    for (size_t i = 0; i < numBins; ++i) {
        normalisedBinPops[i] = (double)generatedBinPopulations[i] / numPoints;
    }

    // Plot a graph of number of points generated and overlay the theoretical thing
    TGraphErrors *testPdf = new TGraphErrors(numBins, binCentres.data(), normalisedBinPops.data(), binWidths.data(), 0);
    testPdf->SetTitle("Test Decay Simulator PDF;time;normalised PDF");

    TGraphErrors *actualPdf = new TGraphErrors(numBins, binCentres.data(), expectedNormalisedBinPop.data(), 0, 0);

    const util::LegendParams_t legendParams = {.x1 = 0.9, .x2 = 0.7, .y1 = 0.7, .y2 = 0.9, .header = ""};

    util::saveObjectsToFile<TGraphErrors>(std::vector<TObject *>{testPdf, actualPdf},
                                          std::vector<std ::string>{"", "SAME"},
                                          std::vector<std::string>{"Expected PDF", "Generated data"},
                                          "testPDF.pdf",
                                          legendParams);

    delete actualPdf;
    delete testPdf;
}

bool SimulatedDecays::isAccepted(const double time, const double uniformVal, bool rightSign)
{
    // a better implementation of this might pass in a function pointer or something
    double funcVal{0};
    if (rightSign) {
        funcVal = _cfRate(time);
    } else {
        funcVal = _dcsRate(time);
    }
    double c = rightSign ? _maxCFRatio : _maxDCSRatio;

    // Check that the RHS of our acc-rej inequality is indeed between 0 and 1
    double rhs = funcVal / (c * _generatingPDF(time));
    if (rhs < 0 || rhs > 1.0) {
        std::string rs = rightSign ? "Right" : "Wrong";
        std::cerr << "For " << rs << " sign decay:" << std::endl;
        std::cerr << "Accept-reject error: f(t)/c*exp(width*t) returned a value " << rhs << " out of range [0, 1]."
                  << std::endl;
        throw D2K3PiException();
    }

    return uniformVal < rhs;
}

void SimulatedDecays::findDcsDecayTimes(size_t numEvents)
{
    // Initialise our vector of wrong-sign decay times to the right length
    WSDecayTimes        = std::vector<double>(numEvents, -1);
    size_t numGenerated = 0;

    while (numGenerated < numEvents) {
        double time = _getRandomTime();
        double y    = _getRandomUniform();

        if (isAccepted(time, y, false)) {
            WSDecayTimes[numGenerated] = time;
            numGenerated++;
        }
    }
}

void SimulatedDecays::findCfDecayTimes(size_t numEvents)
{
    // Initialise our vector of right-sign decay times to the right length
    RSDecayTimes        = std::vector<double>(numEvents, -1);
    size_t numGenerated = 0;

    while (numGenerated < numEvents) {
        double time = _getRandomTime();
        double y    = _getRandomUniform();

        if (isAccepted(time, y, true)) {
            RSDecayTimes[numGenerated] = time;
            numGenerated++;
        }
    }
}

void SimulatedDecays::_setMaxRatios(void)
{
    // Boost's minimising algorithm thing only finds minima, so find the minimum of -1 * rate/generating exponential
    // to find the maximum value.
    auto dcsRatioFunc = [&](double const &x) { return -1 * _dcsRate(x) / _generatingPDF(x); };
    auto cfRatioFunc  = [&](double const &x) { return -1 * _cfRate(x) / _generatingPDF(x); };

    // Find the minima as accurately as possible
    int                       bits   = std::numeric_limits<double>::digits;
    std::pair<double, double> dcsMin = boost::math::tools::brent_find_minima(dcsRatioFunc, _minTime, _maxTime, bits);
    std::pair<double, double> cfMin  = boost::math::tools::brent_find_minima(cfRatioFunc, _minTime, _maxTime, bits);

    // Set _maxDCSRatio to the maximum value of
    _maxDCSRatio = -dcsMin.second;
    _maxCFRatio  = -cfMin.second;
}

void SimulatedDecays::plotRates(const std::vector<double> &timeBinLimits)
{
    // Check tht WSDecayTimes and RSDecayTimes are set
    if (WSDecayTimes.empty()) {
        std::cerr << "WS Decay times not set" << std::endl;
        throw D2K3PiException();
    }
    if (RSDecayTimes.empty()) {
        std::cerr << "RS Decay times not set" << std::endl;
        throw D2K3PiException();
    }

    // Check our time bin limits are sorted
    if (!std::is_sorted(timeBinLimits.begin(), timeBinLimits.end())) {
        std::cerr << "Time bin limits should be sorted" << std::endl;
        throw D2K3PiException();
    }

    size_t numBins = timeBinLimits.size() - 1;
    TH1D * RSHist  = new TH1D("Times", "Right Sign Decay Times;Time/ns;Counts", numBins, timeBinLimits.data());
    TH1D * WSHist  = new TH1D("Times", "Wrong Sign Decay Times;Time/ns;Counts", numBins, timeBinLimits.data());

    RSHist->FillN(RSDecayTimes.size(), RSDecayTimes.data(), nullptr);
    WSHist->FillN(WSDecayTimes.size(), WSDecayTimes.data(), nullptr);
    RSHist->SetStats(false);
    WSHist->SetStats(false);

    util::saveObjectToFile(RSHist, "RSHist.png");
    util::saveObjectToFile(WSHist, "WSHist.png");
    delete RSHist;
    delete WSHist;
}

double SimulatedDecays::_getRandomUniform(void)
{
    return _uniform(*_gen);
}

double SimulatedDecays::maxDCSRatio(void)
{
    return _maxDCSRatio;
}

double SimulatedDecays::maxCFRatio(void)
{
    return _maxCFRatio;
}

#endif // DECAYSIMULATOR_CPP
