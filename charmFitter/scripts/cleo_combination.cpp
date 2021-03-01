#include <algorithm>

#include "CleoCombinationFitter.h"
#include "DecaySimulator.h"
#include "physics.h"

#include "TColor.h"
#include "TEllipse.h"
#include "TGraph2D.h"
#include "TH2D.h"

#include <boost/progress.hpp>

namespace
{

/*
 * Generate N+1 points between -1 and 1
 */
template <size_t N> std::array<double, N + 1> _linspace(void)
{
    std::array<double, N + 1> vals{};

    for (size_t i = 0; i < N + 1; ++i) {
        double val = -1.0 + (i * 2.0) / N;
        vals[i]    = val;
    }

    return vals;
}

/*
 * Generate toy data
 * Returns dcs times, cf times
 */
std::pair<std::vector<double>, std::vector<double>> _toy_data(const FitterUtil::DecayParams_t& params,
                                                              const double                     maxTime)
{
    std::random_device rd;
    std::mt19937       gen{rd()};

    constexpr size_t numCfEvents{2000000};

    SimulatedDecays     DecaySimulator({0.0, maxTime}, params, gen);
    std::vector<double> dcsTimes{DecaySimulator.wsDecayTimes(Phys::numDCSDecays(numCfEvents, params, maxTime))};
    std::vector<double> cfTimes{DecaySimulator.rsDecayTimes(numCfEvents)};

    return {dcsTimes, cfTimes};
}

/*
 * Normalise a vector of chisqs to std devs
 */
std::vector<double> _chisq2sigma(const std::vector<double> chisq)
{
    std::vector<double> sigma = chisq;

    // Subtract off the minimum value
    const double min = *std::min_element(chisq.begin(), chisq.end());
    std::transform(sigma.begin(), sigma.end(), sigma.begin(), [&min](const double x) { return x - min; });

    // Take the sqrt to take chisq -> sigma
    std::transform(sigma.begin(), sigma.end(), sigma.begin(), [](const double x) { return std::sqrt(x); });

    return sigma;
}

/*
 * Format a graph
 */
void _format(TGraph2D& graph, const std::vector<double>& contours)
{
    graph.SetMinimum(0.0);
    graph.SetMaximum(contours.back());
    graph.GetHistogram()->SetContour(contours.size(), contours.data());
    graph.GetZaxis()->SetTitle("\\sigma");
    graph.GetXaxis()->SetTitle("Re(Z)");
    graph.GetYaxis()->SetTitle("Im(Z)");
}

/*
 * Plot a graph
 */
void _plot(TGraph2D& graph, const std::string& path)
{
    // Set colour palette
    constexpr unsigned short numColours{5};
    double                   Length[numColours] = {0.00, 0.50, 0.70, 0.75, 1.0};
    double                   Red[numColours]    = {0.00, 0.00, 1.00, 1.00, 1.00};
    double                   Green[numColours]  = {1.00, 0.00, 0.00, 1.00, 1.00};
    double                   Blue[numColours]   = {0.00, 1.00, 0.00, 1.00, 1.00};
    TColor::CreateGradientColorTable(numColours, Length, Red, Green, Blue, 10000);

    // Set graph properties
    std::vector<double> contours{0.0, 1.0, 2.0, 3.0, 4.0};
    _format(graph, contours);

    // Draw a circle of the allowed Z values
    std::unique_ptr<TEllipse> circle = std::make_unique<TEllipse>(0, 0, 0.56, 0.56);
    circle->SetFillColorAlpha(0, 0);

    // Set canvas properties
    std::unique_ptr<TCanvas> c = std::make_unique<TCanvas>();
    c->cd();
    c->SetWindowSize(1500, 1500);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.15);

    // "Z" options draws a scale bar
    graph.Draw("CONT4");
    graph.SetTitle("");
    circle->Draw();
    c->SaveAs(path.c_str());
}

} // namespace

int main(void)
{
    // Generate allowed values for Z
    constexpr unsigned short N{40};
    const auto               allowedReZ = _linspace<N>();
    const auto               allowedImZ = _linspace<N>();

    // Init empty vectors of reZ, imZ and chisqs
    std::vector<double> reZ{};
    std::vector<double> imZ{};
    std::vector<double> mixingChi2{};
    std::vector<double> cleoChi2{};
    std::vector<double> combinedChi2{};

    // Assign nonsense values of likelihood to 100.0
    constexpr double nonsense{2000.0};

    // Create toy data
    FitterUtil::DecayParams_t params{
        .x     = CharmFitter::WORLD_AVERAGE_X,
        .y     = CharmFitter::WORLD_AVERAGE_Y,
        .r     = 0.055,
        .z_im  = 0.85,
        .z_re  = 0.0,
        .width = 2500.0,
    };
    constexpr double maxTime{0.004};
    const auto [dcsTimes, cfTimes] = _toy_data(params, maxTime);

    // Create fitters
    constexpr unsigned short numBins{25};
    std::vector<double>      binLimits = FitterUtil::exponentialBinLimits(maxTime, params.width, numBins);
    const CLEO::Bin          binNumber{CLEO::Bin::Bin1};
    std::array<double, 6> initialParameterGuess{params.x, params.y, params.r, params.z_im, params.z_re, params.width};
    std::array<double, 6> initialErrorsGuess{1, 1, 1, 1, 1, 1};

    CharmFitter::ConstrainedFitter ConstrainedFitter(binLimits, initialParameterGuess, initialErrorsGuess);
    ConstrainedFitter.fixParameters(std::array<std::string, 3>{"width", "z_im", "z_re"});

    CharmFitter::CLEOCombinationFitter CleoFitter(binLimits, initialParameterGuess, initialErrorsGuess, binNumber);
    CleoFitter.fixParameters(std::array<std::string, 3>{"width", "z_im", "z_re"});

    ConstrainedFitter.addRSPoints(cfTimes, std::vector<double>(cfTimes.size(), 1.0));
    ConstrainedFitter.addWSPoints(dcsTimes, std::vector<double>(dcsTimes.size(), 1.0));

    CleoFitter.addRSPoints(cfTimes, std::vector<double>(cfTimes.size(), 1.0));
    CleoFitter.addWSPoints(dcsTimes, std::vector<double>(dcsTimes.size(), 1.0));

    const auto efficiency = []([[maybe_unused]] const double x) { return 1; };

    // Scan over Z finding the charm fitter values, the CLEO values and the combined values
    boost::progress_display showProgress(allowedReZ.size() * allowedImZ.size());
    for (const auto r : allowedReZ) {
        for (const auto i : allowedImZ) {
            reZ.push_back(r);
            imZ.push_back(i);

            // Construct a set of phsp params
            FitterUtil::DecayParams_t theseParams{params};
            theseParams.z_re = r;
            theseParams.z_im = i;

            // Calculate the CLEO likelihood
            double l = -2 * CLEO::cleoLikelihood(binNumber, theseParams);
            l        = std::isnan(l) ? nonsense : l;
            cleoChi2.push_back(l);

            // Calculate the charm fitter only value
            ConstrainedFitter.setParameter("z_im", i);
            ConstrainedFitter.setParameter("z_re", r);
            mixingChi2.push_back(ConstrainedFitter.fit(efficiency).fitStatistic);

            // Calculate the combined fitter value
            CleoFitter.setParameter("z_im", i);
            CleoFitter.setParameter("z_re", r);
            combinedChi2.push_back(CleoFitter.fit(efficiency).fitStatistic);

            ++showProgress;
        }
    }

    // Normalise each to std devs
    cleoChi2     = _chisq2sigma(cleoChi2);
    mixingChi2   = _chisq2sigma(mixingChi2);
    combinedChi2 = _chisq2sigma(combinedChi2);

    // Make plots of each
    std::unique_ptr<TGraph2D> cleoGraph =
        std::make_unique<TGraph2D>("cleo", "CLEO", cleoChi2.size(), reZ.data(), imZ.data(), cleoChi2.data());
    std::unique_ptr<TGraph2D> mixingGraph =
        std::make_unique<TGraph2D>("mixing", "Mixing", mixingChi2.size(), reZ.data(), imZ.data(), mixingChi2.data());
    std::unique_ptr<TGraph2D> combinedGraph = std::make_unique<TGraph2D>(
        "combined", "Combination", combinedChi2.size(), reZ.data(), imZ.data(), combinedChi2.data());

    _plot(*cleoGraph, "cleo.png");
    _plot(*mixingGraph, "mixing.png");
    _plot(*combinedGraph, "combined.png");

    return 0;
}
