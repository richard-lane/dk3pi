#ifndef MCGENERATOR_CPP
#define MCGENERATOR_CPP

#include <iostream>
#include <random>

#include "../include/D2K3PiError.h"
#include "../include/MCGenerator.h"

MCGenerator::MCGenerator(const double XRangeMin, const double XRangeMax, const double YRangeMin, const double YRangeMax)
{
    _setMinXValue(XRangeMin);
    _setMaxXValue(XRangeMax);

    _setMinYValue(YRangeMin);
    _setMaxYValue(YRangeMax);
}

void MCGenerator::_setMinXValue(const double minValue)
{
    _minX = minValue;
    _setDistributions();
}

void MCGenerator::_setMaxXValue(const double maxValue)
{
    _maxX = maxValue;
    _setDistributions();
}

void MCGenerator::_setMinYValue(const double minValue)
{
    _minY = minValue;
    _setDistributions();
}

void MCGenerator::_setMaxYValue(const double maxValue)
{
    _maxY = maxValue;
    _setDistributions();
}

void MCGenerator::_setDistributions(void)
{
    // Seed our Mersenne twister engine with a random device
    std::random_device rd;
    _gen = std::mt19937(rd());

    // Set distributions to go between our values
    _xDistribution = std::uniform_real_distribution<double>(_minX, _maxX);
    _yDistribution = std::uniform_real_distribution<double>(_minY, _maxY);
}

double MCGenerator::getRandomX(void)
{
    return _xDistribution(_gen);
}

double MCGenerator::getRandomY(void)
{
    return _yDistribution(_gen);
}

bool MCGenerator::isAccepted(const double xVal, const double yVal, const double (*func)(double))
{
    double funcVal = (*func)(xVal);

    // Our Maximum Y value is not large enough to accomodate the function
    if (funcVal > _maxY) {
        std::cout << "Function value " + std::to_string(funcVal) + " is larger than maximum allowed value " +
                         std::to_string(_maxY)
                  << std::endl;

        throw D2K3PiException();
    }

    // Our Maximum X or Y values are smaller than the provided x and y val
    if (xVal > _maxX || xVal < _minX || yVal > _maxY || yVal < _minY) {
        std::cout << "Generated value (" + std::to_string(xVal) + ", " + std::to_string(yVal) +
                         ") is outside of allowed region X(" + std::to_string(_minX) + ", " + std::to_string(_maxX) +
                         "); Y(" + std::to_string(_minY) + ", " + std::to_string(_maxY) + ")"
                  << std::endl;
        throw D2K3PiException();
    }

    return yVal < funcVal;
}

#endif // MCGENERATOR_CPP