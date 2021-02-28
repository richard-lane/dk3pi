/*
 * Make plots showing the generated histograms and fits to them using our different fit parametrisations
 */
#include "ConstrainedFitter.h"
#include "DecaySimulator.h"
#include "PolynomialFitter.h"
#include "UnconstrainedFitter.h"
#include "fitterUtil.h"
#include "physics.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1D.h>

namespace
{
}

int main()
{
    // Options and parameters
    const FitterUtil::DecayParams_t decayParams{
        .x     = CharmFitter::WORLD_AVERAGE_X,
        .y     = CharmFitter::WORLD_AVERAGE_Y,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = 2439,
    };
    constexpr double maxTime{0.004};
    constexpr size_t numRSEvents{10000000};
    const size_t     numWSEvents = static_cast<size_t>(Phys::numDCSDecays(numRSEvents, decayParams, maxTime));

    const auto efficiency = [](const double) { return 1.0; };

    constexpr size_t          numBins{50};
    const std::vector<double> expBinLimits{FitterUtil::exponentialBinLimits(maxTime, decayParams.width, numBins)};
    std::vector<double>       linearBinLimits(numBins + 1);
    for (size_t i = 0; i <= numBins; ++i) {
        linearBinLimits[i] = maxTime * i / numBins;
    }

    // Generate some decay times
    std::random_device  rd;
    std::mt19937        gen(rd());
    SimulatedDecays     DecaySimulator({0.0, maxTime}, decayParams, gen);
    std::vector<double> rsTimes = DecaySimulator.rsDecayTimes(numRSEvents);
    std::vector<double> wsTimes = DecaySimulator.wsDecayTimes(numWSEvents);

    // Histogram them
    TH1D rsHist("RS", "RS", numBins, linearBinLimits.data());
    TH1D wsHist("WS", "WS", numBins, linearBinLimits.data());
    rsHist.FillN(numRSEvents, rsTimes.data(), nullptr);
    wsHist.FillN(numWSEvents, wsTimes.data(), nullptr);
    rsHist.SetStats(false);
    wsHist.SetStats(false);
    wsHist.SetLineColor(kRed);

    // Perform fits to them
    CharmFitter::CharmPolynomialFitter PolyFitter(expBinLimits, {1, 1, 1}, {1, 1, 1}, decayParams.width);
    PolyFitter.addRSPoints(rsTimes, std::vector(numRSEvents, 1.0));
    PolyFitter.addWSPoints(wsTimes, std::vector(numWSEvents, 1.0));
    auto PolyFunction{PolyFitter.fit(efficiency).bestFitFunction};

    CharmFitter::ConstrainedFitter ConstrainedFitter(
        expBinLimits,
        {decayParams.x, decayParams.y, decayParams.r, decayParams.z_im, decayParams.z_re, decayParams.width},
        {1, 1, 1, 1, 1, 1});
    ConstrainedFitter.addRSPoints(rsTimes, std::vector(numRSEvents, 1.0));
    ConstrainedFitter.addWSPoints(wsTimes, std::vector(numWSEvents, 1.0));
    ConstrainedFitter.fixParameters(std::array{"width", "z_re", "z_im"});
    auto ConstrainedFunction{ConstrainedFitter.fit(efficiency).bestFitFunction};
    ConstrainedFunction.SetLineColor(kBlue);

    // Plot the fits
    TCanvas Canvas("", "", 1200, 1600);
    Canvas.Divide(1, 3, 0.0, 0.0);
    Canvas.cd(1);
    const auto ratios{PolyFitter.ratios()};
    const auto errors{PolyFitter.errors()};
    const auto binCentres{PolyFitter.getBinCentres()};
    auto       binErrors{PolyFitter.getBinWidths()};
    std::transform(binErrors.begin(), binErrors.end(), binErrors.begin(), [](const double x) { return x / 2.0; });
    TGraphErrors data(numBins, binCentres.data(), ratios.data(), binErrors.data(), errors.data());

    data.Draw("AP");
    PolyFunction.Draw("SAME");
    ConstrainedFunction.Draw("SAME");

    Canvas.cd(2);
    rsHist.Draw();

    Canvas.cd(3);
    wsHist.Draw();

    Canvas.SaveAs("plot.png");

    return 0;
}
