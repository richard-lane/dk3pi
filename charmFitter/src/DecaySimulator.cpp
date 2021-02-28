#include "DecaySimulator.h"
#include "fitterUtil.h"
#include "physics.h"

SimulatedDecays::SimulatedDecays(const std::pair<double, double>& timeDomain,
                                 const FitterUtil::DecayParams_t& decayParams,
                                 std::mt19937&                    rng)
    : _minTime(timeDomain.first), _maxTime(timeDomain.second), _decayParams(decayParams), _gen(&rng)
{
    _uniform    = std::uniform_real_distribution<double>(0.0, 1.0);
    _maxWSRatio = _findMaxRatio();
}

double SimulatedDecays::rsPoint()
{
    double x = _uniform(*_gen);
    double z = 1 - std::exp(-1 * _decayParams.width * _maxTime);
    return (-1 / _decayParams.width) * std::log(1 - z * x);
}

double SimulatedDecays::wsPoint()
{
    // Loop until we accept a point
    while (true) {
        // Generate a time from our exponential distribution
        double exponentialTime{this->rsPoint()};

        // Check whether it is accepted
        if (_isAccepted(exponentialTime)) {
            return exponentialTime;
        }
    }
}

std::vector<double> SimulatedDecays::rsDecayTimes(const size_t numPoints)
{
    std::vector<double> times(numPoints);

    for (size_t i = 0; i < numPoints; ++i) {
        times[i] = rsPoint();
    }

    return times;
}

std::vector<double> SimulatedDecays::wsDecayTimes(const size_t numPoints)
{
    std::vector<double> times(numPoints);

    for (size_t i = 0; i < numPoints; ++i) {
        times[i] = wsPoint();
    }

    return times;
}

double SimulatedDecays::_rsRate(const double time) const
{
    return std::exp(-_decayParams.width * time);
}

double SimulatedDecays::_wsRate(const double time) const
{
    return Phys::rateRatio(time, _decayParams) * _rsRate(time);
}

double SimulatedDecays::_findMaxRatio() const
{
    // Find whether our ratio is higher at mintime or maxtime
    const double minTimeRatio{Phys::rateRatio(_minTime, _decayParams)};
    const double maxTimeRatio{Phys::rateRatio(_maxTime, _decayParams)};
    double       maxRatio = minTimeRatio > maxTimeRatio ? minTimeRatio : maxTimeRatio;

    // Find whether the potential turning point of (a + bt + ct^2) is in the range of allowed times
    // If it is, find whether the ratio there is larger than at _minTime or _maxTime
    const std::array abc{Phys::expectedParams(_decayParams)};
    const double     tpTime{-abc[1] / (2 * abc[2])};
    if (_minTime <= tpTime && tpTime <= _maxTime) {
        maxRatio = Phys::rateRatio(tpTime, _decayParams) > maxRatio ? Phys::rateRatio(tpTime, _decayParams) : maxRatio;
    }

    return maxRatio;
}

bool SimulatedDecays::_isAccepted(const double time)
{
    // Find a random number uniformly between 0 and 1, then accept it if it is <= (a + bt + ct^2)/maxRatio
    return _uniform(*_gen) < Phys::rateRatio(time, _decayParams) / _maxWSRatio;
}
