#include <boost/test/unit_test.hpp>

#include "../pull_study/DecaySimulator.h"

BOOST_AUTO_TEST_SUITE(simulator_IT)

/*
 * Integration-style test that the simulator correctly generates x^2 and x
 */
BOOST_AUTO_TEST_CASE(it_simulator)
{
    double              maxTime     = 1.;
    size_t              numTimeBins = 50;
    std::vector<double> timeBinLimits{};
    for (size_t i = 0; i < numTimeBins + 1; ++i) {
        timeBinLimits.push_back(i * (maxTime / numTimeBins));
    }

    // Lambdas for our functions to generate
    auto linear    = [](double x) { return x; };     // CF
    auto quadratic = [](double x) { return x * x; }; // DCS

    // PDF and lambda to generate it
    std::random_device                     rd;
    std::shared_ptr<std::mt19937>          _gen = std::make_shared<std::mt19937>(rd());
    std::uniform_real_distribution<double> uniform;

    // Use e^-x as our PDF
    auto pdf           = [&](double x) { return std::exp(-1 * x) / (1 - std::exp(-1 * maxTime)); };
    auto numberFromPDF = [&](void) {
        double x = uniform(*_gen);
        double z = 1 - std::exp(-1 * maxTime);
        return -1 * std::log(1 - z * x);
    };

    SimulatedDecays MyDecays =
        SimulatedDecays(numberFromPDF, pdf, linear, quadratic, std::make_pair(0.0, maxTime), _gen);

    size_t N = 1e6;
    MyDecays.findCfDecayTimes(N);
    MyDecays.findDcsDecayTimes(N);

    // Bin decays
    std::vector<size_t> numCF  = util::binVector(MyDecays.RSDecayTimes, timeBinLimits);
    std::vector<size_t> numDCS = util::binVector(MyDecays.WSDecayTimes, timeBinLimits);

    // Find expected numbers of decays in each bin
    std::vector<double> expectedNumCF(numTimeBins, -1);
    std::vector<double> expectedNumDCS(numTimeBins, -1);
    for (size_t i = 0; i < numTimeBins; ++i) {
        expectedNumCF[i] = N * 0.5 *
                           (timeBinLimits[i + 1] * timeBinLimits[i + 1] - timeBinLimits[i] * timeBinLimits[i]) /
                           (0.5 * maxTime * maxTime);
        expectedNumDCS[i] = N * (1. / 3.) * (std::pow(timeBinLimits[i + 1], 3) - std::pow(timeBinLimits[i], 3)) /
                            ((1. / 3.) * maxTime * maxTime * maxTime);
    }

    // Find chi squared values
    double cfChiSq{0};
    double dcsChiSq{0};

    for (size_t i = 0; i < numTimeBins; ++i) {
        cfChiSq += std::pow(expectedNumCF[i] - (double)numCF[i], 2) / expectedNumCF[i];
        dcsChiSq += std::pow(expectedNumDCS[i] - (double)numDCS[i], 2) / expectedNumDCS[i];
    }

    // Expect chi squared to be distributed with mean of 50; if the model is wrong expect huge deviation
    BOOST_CHECK(cfChiSq < 100 && cfChiSq > 1);
    BOOST_CHECK(dcsChiSq < 100 && dcsChiSq > 1);
    std::cout << "chi squares: " << cfChiSq << " " << dcsChiSq << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()