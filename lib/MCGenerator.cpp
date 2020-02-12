#ifndef MCGENERATOR_CPP
#define MCGENERATOR_CPP

#include <iostream>
#include <random>

#include "D2K3PiError.h"
#include "MCGenerator.h"

MCGenerator::MCGenerator(const std::pair<double, double> &xRange, const std::pair<double, double> &yRange)
{
    setXRange(xRange);
    setYRange(yRange);
}

void MCGenerator::setXRange(const std::pair<double, double> &xRange)
{
    if (xRange.first > xRange.second) {
        std::cerr << "Provided X values (" + std::to_string(xRange.first) + ", " + std::to_string(xRange.second) +
                         ") do not define a range. Have you created your std::pair backwards?"
                  << std::endl;
        throw D2K3PiException();
    }
    _minX = xRange.first;
    _maxX = xRange.second;
    setDistributions();
}

void MCGenerator::setYRange(const std::pair<double, double> &yRange)
{
    if (yRange.first > yRange.second) {
        std::cerr << "Provided Y values (" + std::to_string(yRange.first) + ", " + std::to_string(yRange.second) +
                         ") do not define a range. Have you created your std::pair backwards?"
                  << std::endl;
        throw D2K3PiException();
    }
    _minY = yRange.first;
    _maxY = yRange.second;
    setDistributions();
}

void MCGenerator::setDistributions(void)
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

bool MCGenerator::isAccepted(const double xVal, const double yVal, double (*func)(double))
{
    double funcVal = (*func)(xVal);

    // Our Maximum Y value is not large enough to accomodate the function
    if (funcVal > _maxY) {
        std::cerr<< "Function value " + std::to_string(funcVal) + " is larger than maximum allowed value " +
                         std::to_string(_maxY)
                  << std::endl;

        throw D2K3PiException();
    }

    // Our Maximum X or Y values are smaller than the provided x and y val
    if (xVal > _maxX || xVal < _minX || yVal > _maxY || yVal < _minY) {
        std::cerr<< "Generated value (" + std::to_string(xVal) + ", " + std::to_string(yVal) +
                         ") is outside of allowed region X(" + std::to_string(_minX) + ", " + std::to_string(_maxX) +
                         "); Y(" + std::to_string(_minY) + ", " + std::to_string(_maxY) + ")"
                  << std::endl;
        throw D2K3PiException();
    }

    return yVal < funcVal;
}

#endif // MCGENERATOR_CPP