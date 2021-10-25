#include <random>
#include <iostream>

#include "simulator.h"
#include "fitterUtil.h"
#include "physics.h"
#include "DecaySimulator.h"

static FitterUtil::DecayParams_t _arr2params(const std::array<double, 6>& decayParamsArr) {
    return FitterUtil::DecayParams_t{decayParamsArr[0],
                                     decayParamsArr[1],
                                     decayParamsArr[2],
                                     decayParamsArr[3],
                                     decayParamsArr[4],
                                     decayParamsArr[5]};
}

std::array<double, 3> expectedParamsBinding(const std::array<double, 6>& decayParamsArr) {
    const auto decayParams = _arr2params(decayParamsArr);
    return Phys::expectedParams(decayParams);
}

std::pair<std::vector<double>, std::vector<double>>
simulate(const size_t n,
         const std::array<double, 6>& decayParamsArr,
         const double maxTime,
         const uint32_t seed) {
    // Convert array to decay params
    const auto decayParams = _arr2params(decayParamsArr);

    // Find how many WS evts to make
    double numWS = Phys::numDCSDecays(n, decayParams, maxTime);
    //std::cout << "Creating " << n << " RS and " << numWS << "WS evts" << std::endl;

    // Create simulator
    std::mt19937    rng{seed};
    SimulatedDecays Simulator({0.0, maxTime}, decayParams, rng);

    // Simulate evts
    //std::cout << "creating RS..." << std::endl;
    std::vector<double> rsTimes = Simulator.rsDecayTimes(n);

    //std::cout << "creating WS..." << std::endl;
    std::vector<double> wsTimes = Simulator.wsDecayTimes(numWS);

    // Return
    return {rsTimes, wsTimes};
}

