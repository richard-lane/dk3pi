#include "DecaySimulator.h"
#include "PolynomialFitter.h"
#include "physics.h"
#include "util.h"

#include <TGraphErrors.h>
#include <TH1D.h>

#include <boost/progress.hpp>

namespace
{
/*
 * Find limits describing N bins
 */
template <size_t N> std::array<double, N + 1> binLimits(const double maxTime)
{
    std::array<double, N + 1> bins;
    for (size_t i = 0; i <= N; ++i) {
        bins[i] = i * maxTime / N;
    }
    return bins;
}

/*
 * Create histograms of the times generated
 */
void plot_hists()
{
    // Some params
    constexpr size_t numBins{25};
    constexpr double maxTime{0.004};
    size_t           numEvents{10000000};

    FitterUtil::DecayParams_t decayParams{
        .x = 0.0039, .y = 0.0065, .r = 0.055, .z_im = -0.2956, .z_re = 0.7609, .width = 2400};

    // Create a RNG
    std::random_device rd{};
    std::mt19937       gen(rd());

    // Create ROOT histograms
    auto timeBins{binLimits<numBins>(maxTime)};
    TH1D rsHist("RS Simulation", "RS Simulation", numBins, timeBins.data());
    TH1D wsHist("WS Simulation", "WS Simulation", numBins, timeBins.data());

    // Create a decay simulator
    SimulatedDecays Simulator({0.0, maxTime}, decayParams, gen);

    // Progress bar
    boost::progress_display progressBar{numEvents};

    // Fill these histograms
    for (size_t i = 0; i < numEvents; ++i) {
        rsHist.Fill(Simulator.rsPoint());
        wsHist.Fill(Simulator.wsPoint());
        ++progressBar;
    }
    wsHist.SetLineColor(kRed);

    // Plot them on the same canvas, and on different canvases
    util::LegendParams_t legend{.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = ""};
    util::saveObjectsToFile<TH1D>({&rsHist, &wsHist}, {"", "SAME"}, {"RS", "WS"}, "sim.png", legend);
}

/*
 * If we produce equal numbers of WS/RS events, we expect the ratio to look like
 *     k(a + bt + ct^2)
 * This function finds the coefficienct k
 */
inline double wsRSRatio(const FitterUtil::DecayParams_t& decayParams, const double maxTime)
{
    // Create some aliases to make the expression more readable
    const auto    abc{Phys::expectedParams(decayParams)};
    const double& a{abc[0]};
    const double& b{abc[1]};
    const double& c{abc[2]};
    const double& w{decayParams.width};
    const double& T{maxTime};

    double denominator =
        w * (a * w + b) + 2 * c - std::exp(-w * T) * (w * (a * w + b * w * T + b) + c * (w * T * (w * T + 2) + 2));
    double numerator = w * w * (1 - std::exp(-w * T));

    return numerator / denominator;
}

/*
 * Find the expected ratio WS/RS (if equal numbers of events generated) at a given time
 */
inline double expectedRatio(const FitterUtil::DecayParams_t& decayParams, const double maxTime, const double time)
{
    const auto abc{Phys::expectedParams(decayParams)};
    return wsRSRatio(decayParams, maxTime) * (abc[0] + abc[1] * time + abc[2] * time * time);
}

/*
 * Replace non-finite elements of a collection with 0, in place
 */
template <typename T> void removeNonFinite(T& c)
{
    for (auto& x : c) {
        if (!std::isfinite(x)) {
            x = 0;
        }
    }
}

/*
 * Check that the ratio of generated decays follows the right shape
 */
void ratio_check()
{
    // Some params
    // Need narrow bins so we don't lose out by integrating over them to find expected events
    constexpr size_t numBins{250};
    constexpr double maxTime{0.008};
    size_t           numEvents{10000000};

    FitterUtil::DecayParams_t decayParams{
        .x = 0.0039, .y = 0.0065, .r = 0.055, .z_im = -0.2956, .z_re = 0.7609, .width = 2400};

    // Create a RNG
    std::random_device rd{};
    std::mt19937       gen(rd());

    // Create a decay simulator
    SimulatedDecays Simulator({0.0, maxTime}, decayParams, gen);

    // Create a Fitter that we'll just use to do the time binning
    auto                               timeBins{binLimits<numBins>(maxTime)};
    CharmFitter::CharmPolynomialFitter Fitter({timeBins.begin(), timeBins.end()}, {1, 1, 1}, {1, 1, 1});

    // Progress bar
    boost::progress_display progressBar{numEvents};

    // Generate points + add them to the histograms stored in the fitter
    for (size_t i = 0; i < numEvents; ++i) {
        Fitter.addRSPoint(Simulator.rsPoint());
        Fitter.addWSPoint(Simulator.wsPoint());
        ++progressBar;
    }

    // Find what we expect the ratio at each bin centre to be
    const std::vector<double> binCentres{Fitter.getBinCentres()};
    std::vector<double>       expectedRatios{binCentres};
    std::transform(binCentres.begin(), binCentres.end(), expectedRatios.begin(), [&decayParams, &maxTime](double t) {
        return expectedRatio(decayParams, maxTime, t);
    });

    // Create ROOT graphs
    std::vector<double> ratios{Fitter.ratios()};
    std::vector<double> errors{Fitter.errors()};
    removeNonFinite(ratios);
    removeNonFinite(errors);
    TGraphErrors graph(numBins, binCentres.data(), ratios.data(), nullptr, errors.data());
    TGraphErrors expected(numBins, binCentres.data(), expectedRatios.data(), nullptr, nullptr);
    expected.SetLineColor(kRed);

    util::LegendParams_t legend{.x1 = 0.9, .x2 = 0.7, .y1 = 0.1, .y2 = 0.3, .header = ""};
    util::saveObjectsToFile<TGraph>(
        {&graph, &expected}, {"AP*", "L SAME"}, {"Generated", "Expected"}, "ratios.png", legend);
}

} // namespace

int main()
{
    plot_hists();
    ratio_check();
}
