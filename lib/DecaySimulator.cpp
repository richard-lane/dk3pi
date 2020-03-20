#ifndef DECAYSIMULATOR_CPP
#define DECAYSIMULATOR_CPP

#include <cmath>
#include <iostream>
#include <utility>

#include <TH1D.h>

#include "D2K3PiError.h"
#include "DecaySimulator.h"
#include "util.h"

SimulatedDecays::SimulatedDecays(const std::pair<double, double> &timeRange,
                                 const std::pair<double, double> &decayRateRange,
                                 const DecayParams_t &            DecayParams)
{
    // The parent class implementations will suffice to set the allowed ranges of values.
    setXRange(timeRange);
    setYRange(decayRateRange);

    _DecayParams = DecayParams;
}

bool SimulatedDecays::isAccepted(const double xVal, const double yVal, bool rightSign)
{
    // Note: a better implementation of this might use the parent class implementation of isAccepted,
    // but I didn't want to do that because of issues with pointers to member functions being different
    // from normal function pointers :(
    double funcVal{0};
    if (rightSign) {
        funcVal = _rightSignDecayRate(xVal);
    } else {
        funcVal = _wrongSignDecayRate(xVal);
    }

    // Our Maximum Y value is not large enough to accomodate the function
    if (funcVal > _maxY) {
        std::cerr << "Function value " + std::to_string(funcVal) + " is larger than maximum allowed value " +
                         std::to_string(_maxY)
                  << std::endl;

        throw D2K3PiException();
    }

    // Our Maximum X or Y values are smaller than the provided x and y val
    if (xVal > _maxX || xVal < _minX || yVal > _maxY || yVal < _minY) {
        std::cerr << "Generated value (" + std::to_string(xVal) + ", " + std::to_string(yVal) +
                         ") is outside of allowed region X(" + std::to_string(_minX) + ", " + std::to_string(_maxX) +
                         "); Y(" + std::to_string(_minY) + ", " + std::to_string(_maxY) + ")"
                  << std::endl;
        throw D2K3PiException();
    }

    return yVal < funcVal;
}

void SimulatedDecays::findDcsDecayTimes(size_t numEvents)
{
    // Initialise our vector of wrong-sign decay times to the right length
    WSDecayTimes        = std::vector<double>(numEvents, -1);
    size_t numGenerated = 0;

    while (numGenerated < numEvents) {
        double time  = getRandomX();
        double ratio = getRandomY();

        if (isAccepted(time, ratio, false)) {
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
        double time  = getRandomX();
        double ratio = getRandomY();

        if (isAccepted(time, ratio, true)) {
            RSDecayTimes[numGenerated] = time;
            numGenerated++;
        }
    }
}

double SimulatedDecays::_rightSignDecayRate(const double time)
{
    return exp(-1.0 * _DecayParams.width * time);
}

double SimulatedDecays::_wrongSignDecayRate(const double time)
{
    // Write the decay rate as (a + bt + ct^2)e^(-gamma*t) (ignoring overall factor of B^2 that has been taken out)
    double a = pow(_DecayParams.r, 2);
    double b =
        _DecayParams.r * (_DecayParams.y * _DecayParams.z_re + _DecayParams.x * _DecayParams.z_im) * _DecayParams.width;
    double c = 0.25 * (pow(_DecayParams.x, 2) + pow(_DecayParams.y, 2)) * pow(_DecayParams.width, 2);

    return (a + b * time + c * pow(time, 2)) * exp(-1.0 * _DecayParams.width * time);
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

    // Check our time bin limits are sorted and cover the entire range of time values
    size_t numBins = timeBinLimits.size() - 1;
    if (!std::is_sorted(timeBinLimits.begin(), timeBinLimits.end())) {
        std::cerr << "Time bin limits should be sorted" << std::endl;
        throw D2K3PiException();
    }
    if (timeBinLimits[0] > _minX || timeBinLimits[numBins] < _maxX) {
        std::cerr << "Time bin limits do not cover entire range of possible time values: " << _minX << ", " << _maxX
                  << std::endl;
        throw D2K3PiException();
    }

    TH1D *RSHist = new TH1D("Test accept-reject, RS", "", numBins, timeBinLimits.data());
    TH1D *WSHist = new TH1D("Test accept-reject, WS", "", numBins, timeBinLimits.data());

    for (auto it = RSDecayTimes.begin(); it != RSDecayTimes.end(); ++it) {
        RSHist->Fill(*it);
    }

    for (auto it = WSDecayTimes.begin(); it != WSDecayTimes.end(); ++it) {
        WSHist->Fill(*it);
    }

    util::saveObjectToFile(RSHist, "RSHist.pdf");
    util::saveObjectToFile(WSHist, "WSHist.pdf");
    delete RSHist;
    delete WSHist;
}

#endif // DECAYSIMULATOR_CPP
