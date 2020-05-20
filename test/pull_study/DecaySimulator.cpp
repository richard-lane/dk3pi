#ifndef DECAYSIMULATOR_CPP
#define DECAYSIMULATOR_CPP

#include <cmath>
#include <iostream>
#include <utility>

#include <TH1D.h>

#include "D2K3PiError.h"
#include "DecaySimulator.h"
#include "PullStudyHelpers.h"
#include "physics.h"
#include "util.h"

SimulatedDecays::SimulatedDecays(const double maxTime, const DecayParams_t &DecayParams) : _maxTime(maxTime)
{
    _DecayParams = DecayParams;
    _setMaxDCSRatio();

    // Set our generator to a new random device
    std::random_device rd;
    _gen = std::mt19937(rd());
}

bool SimulatedDecays::isAccepted(const double time, const double uniformVal, bool rightSign)
{
    // At the moment the efficiency function is a step function at tau
    double tau = 1 / _DecayParams.width;
    // Note: a better implementation of this might pass in a function pointer or something
    double funcVal{0};
    if (rightSign) {
        funcVal = Phys::cfRateWithEfficiency(time, _DecayParams, tau);
    } else {
        funcVal = Phys::dcsRateWithEfficiency(time, _DecayParams, tau);
    }
    double c = rightSign ? 1.0 : _maxDCSRatio;

    // Check that the RHS of our acc-rej inequality is indeed between 0 and 1
    double rhs = funcVal / (c * std::exp(-1.0 * _DecayParams.width * time));
    if (rhs < 0 || rhs > 1.0) {
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

void SimulatedDecays::_setMaxDCSRatio(void)
{
    std::vector<double> params = util::expectedParams(_DecayParams);
    double              a      = params[0];
    double              b      = params[1];
    double              c      = params[2];

    _maxDCSRatio = a > a + b * _maxTime + c * _maxTime * _maxTime ? a : a + b * _maxTime + c * _maxTime * _maxTime;
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

    util::saveObjectToFile(RSHist, "RSHist.png");
    util::saveObjectToFile(WSHist, "WSHist.png");
    delete RSHist;
    delete WSHist;
}

double SimulatedDecays::_getRandomUniform(void)
{
    return _uniform(_gen);
}

double SimulatedDecays::_getRandomTime(void)
{
    // Get a random number from our uniform distribution and use an analytical formula to convert it to one from an
    // exponential distribution up to a maximum time.
    double x = _getRandomUniform();
    double z = 1 - std::exp(-1 * _DecayParams.width * _maxTime);
    return (-1 / _DecayParams.width) * std::log(1 - z * x);
}

#endif // DECAYSIMULATOR_CPP
