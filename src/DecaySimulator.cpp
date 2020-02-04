#ifndef DECAYSIMULATOR_CPP
#define DECAYSIMULATOR_CPP

#include <cmath>
#include <iostream>
#include <utility>

#include <TGraph.h>

#include "../include/D2K3PiError.h"
#include "../include/DecaySimulator.h"
#include "../include/util.h"

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
    double c = 0.25 * pow(_DecayParams.x, 2) * pow(_DecayParams.y, 2) * pow(_DecayParams.width, 2);

    return (a + b * time + c * pow(time, 2)) * exp(-1.0 * _DecayParams.width * time);
}

void SimulatedDecays::plotRates(void)
{
    size_t N = 1000;
    // Create vectors of times and both rates
    std::vector<double> times  = std::vector<double>(N, 0);
    std::vector<double> rsRate = std::vector<double>(N, 0);
    std::vector<double> wsRate = std::vector<double>(N, 0);

    for (size_t i = 0; i < N; ++i) {
        double time = 0.000002 * i;
        times[i]    = time;
        rsRate[i]   = _rightSignDecayRate(time);
        wsRate[i]   = _wrongSignDecayRate(time);
    }

    TGraph *rsGraph = new TGraph(N, times.data(), rsRate.data());
    TGraph *wsGraph = new TGraph(N, times.data(), wsRate.data());
    util::saveToFile(rsGraph, "rs.pdf", "AP");
    util::saveToFile(wsGraph, "ws.pdf", "AP");

    delete rsGraph;
    delete wsGraph;
}

#endif // DECAYSIMULATOR_CPP
