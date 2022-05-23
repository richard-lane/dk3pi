#include <algorithm>

#include "BesCombinationFitter.h"
#include "DecaySimulator.h"
#include "physics.h"

#include "TColor.h"
#include "TEllipse.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TMarker.h"

#include <boost/progress.hpp>

namespace
{

/*
 * Generate N points between -x and x, inclusive
 */
template <size_t N, typename T> std::array<double, N> _linspace(const T x)
{
    std::array<double, N> vals{};

    for (size_t i = 0; i < N; ++i) {
        double val = x * (-1.0 + (i * 2.0 / (N - 1)));
        vals[i]    = val;
    }

    return vals;
}

/*
 * Generate toy data
 * Returns dcs times, cf times
 */
std::pair<std::vector<double>, std::vector<double>>
_toy_data(const FitterUtil::DecayParams_t& params, const double maxTime, std::mt19937& gen)
{
    constexpr size_t numCfEvents{2000000};
    const size_t     numDcsEvents{static_cast<size_t>(Phys::numDCSDecays(numCfEvents, params, maxTime))};

    std::cout << numCfEvents << " CF events\t" << numDcsEvents << " DCS events" << std::endl;

    SimulatedDecays     DecaySimulator({0.0, maxTime}, params, gen);
    std::vector<double> dcsTimes{DecaySimulator.wsDecayTimes(numDcsEvents)};
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
 * Plot a graph, and a star denoting a point
 */
void _plot(TGraph2D& graph, const std::string& path, const std::pair<double, double>& point)
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

    // Set canvas properties
    std::unique_ptr<TCanvas> c = std::make_unique<TCanvas>();
    c->cd();
    c->SetWindowSize(1500, 1500);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.15);

    // Create a pad to draw stuff on
    std::unique_ptr<TPad> pad = std::make_unique<TPad>("nullPad", "", 0, 0, 1, 1);
    pad->Draw();
    pad->cd();

    graph.Draw("CONT4"); // "Z" options draws a scale bar
    graph.SetTitle("");

    // Plot our point
    std::unique_ptr<TMarker> pointPlot = std::make_unique<TMarker>(point.first, point.second, 0);
    pointPlot->SetMarkerStyle(29);
    pointPlot->SetMarkerSize(3);
    pointPlot->SetMarkerColor(kYellow);

    // The way to do this is to create a transparent pad on top of our current pad
    std::unique_ptr<TPad> nullPad = std::make_unique<TPad>("nullPad", "", 0, 0, 1, 1);
    nullPad->SetFillStyle(0);
    nullPad->SetFrameFillStyle(0);
    nullPad->Draw();
    nullPad->cd();

    double bm = pad->GetBottomMargin();
    double lm = pad->GetLeftMargin();
    double rm = pad->GetRightMargin();
    double to = pad->GetTopMargin();
    double x1 = graph.GetXaxis()->GetXmin();
    double yf = graph.GetYaxis()->GetXmin();
    double x2 = graph.GetXaxis()->GetXmax();
    double y2 = graph.GetYaxis()->GetXmax();

    double Xa = (x2 - x1) / (1 - lm - rm) - (x2 - x1);
    double Ya = (y2 - yf) / (1 - bm - to) - (y2 - yf);
    double LM = Xa * (lm / (lm + rm));
    double RM = Xa * (rm / (lm + rm));
    double BM = Ya * (bm / (bm + to));
    double TM = Ya * (to / (bm + to));

    nullPad->Range(x1 - LM, yf - BM, x2 + RM, y2 + TM);

    pointPlot->Draw();

    // Draw a circle of the allowed Z values
    std::unique_ptr<TEllipse> circle = std::make_unique<TEllipse>(0, 0, 1.0, 1.0);
    circle->SetFillColorAlpha(0, 0);
    circle->Draw();

    c->SaveAs(path.c_str());
}

/*
 * Appends to vectors in place
 */
template <size_t N>
void _besScan(const BES::Bin                   binNumber,
              const FitterUtil::DecayParams_t& params,
              const double                     nonsense,
              const std::array<double, N>&     allowedReZ,
              const std::array<double, N>&     allowedImZ,
              std::vector<double>&             reZ,
              std::vector<double>&             imZ,
              std::vector<double>&             besChi2)
{
    for (const auto r : allowedReZ) {
        for (const auto i : allowedImZ) {
            reZ.push_back(r);
            imZ.push_back(i);

            // Construct a set of phsp params
            FitterUtil::DecayParams_t theseParams{params};
            theseParams.z_re = r;
            theseParams.z_im = i;

            // Calculate the BES chi2
            double c = BES::besLikelihood(binNumber, theseParams);
            c        = std::isnan(c) ? nonsense : c;
            besChi2.push_back(c);
        }
    }
}

template <size_t N>
void _makePlots(const BES::Bin                   binNumber,
                const double                     maxTime,
                const FitterUtil::DecayParams_t& params,
                const double                     nonsense,
                const std::array<double, N>&     allowedReZ,
                const std::array<double, N>&     allowedImZ,
                const std::string&               prefix)
{
    auto paramsCopy = params;

    // Init empty vectors of reZ, imZ and chisqs
    std::vector<double> reZ{};
    std::vector<double> imZ{};
    std::vector<double> mixingChi2{};
    std::vector<double> besChi2{};
    std::vector<double> combinedChi2{};

    // Scan over Z finding the BES chi2 values
    _besScan(binNumber, paramsCopy, nonsense, allowedReZ, allowedImZ, reZ, imZ, besChi2);

    // Find where the minimum BES chi2 is
    unsigned minBesChi2Index = std::min_element(besChi2.begin(), besChi2.end()) - besChi2.begin();

    // Create toy data with Z at this minimum
    std::random_device rd;
    std::mt19937       gen{rd()};

    paramsCopy.z_im                = imZ[minBesChi2Index];
    paramsCopy.z_re                = reZ[minBesChi2Index];
    const auto [dcsTimes, cfTimes] = _toy_data(paramsCopy, maxTime, gen);

    // Create fitters
    constexpr unsigned short numBins{25};
    std::vector<double>      binLimits = FitterUtil::exponentialBinLimits(maxTime, paramsCopy.width, numBins);
    std::array<double, 6>    initialParameterGuess{
        paramsCopy.x, paramsCopy.y, paramsCopy.r, paramsCopy.z_im, paramsCopy.z_re, paramsCopy.width};
    std::array<double, 6> initialErrorsGuess{1, 1, 1, 1, 1, 1};

    CharmFitter::ConstrainedFitter ConstrainedFitter(binLimits, initialParameterGuess, initialErrorsGuess);
    ConstrainedFitter.fixParameters(std::array<std::string, 3>{"width", "z_im", "z_re"});

    CharmFitter::BESCombinationFitter BesFitter(binLimits, initialParameterGuess, initialErrorsGuess, binNumber);
    BesFitter.fixParameters(std::array<std::string, 3>{"width", "z_im", "z_re"});

    ConstrainedFitter.addRSPoints(cfTimes, std::vector<double>(cfTimes.size(), 1.0));
    ConstrainedFitter.addWSPoints(dcsTimes, std::vector<double>(dcsTimes.size(), 1.0));

    BesFitter.addRSPoints(cfTimes, std::vector<double>(cfTimes.size(), 1.0));
    BesFitter.addWSPoints(dcsTimes, std::vector<double>(dcsTimes.size(), 1.0));

    const auto efficiency = []([[maybe_unused]] const double x) { return 1; };

    // Scan over Z finding the charm fitter values, the BES values and the combined values
    boost::progress_display showProgress(allowedReZ.size() * allowedImZ.size());
    for (const auto r : allowedReZ) {
        for (const auto i : allowedImZ) {
            // Calculate the charm fitter only value
            ConstrainedFitter.setParameter("z_im", i);
            ConstrainedFitter.setParameter("z_re", r);
            mixingChi2.push_back(ConstrainedFitter.fit(efficiency).fitStatistic);

            // Calculate the combined fitter value
            BesFitter.setParameter("z_im", i);
            BesFitter.setParameter("z_re", r);
            combinedChi2.push_back(BesFitter.fit(efficiency).fitStatistic);

            ++showProgress;
        }
    }

    // Normalise each to std devs
    std::cout << "Min BES chi2" << reZ[minBesChi2Index] << " " << imZ[minBesChi2Index] << std::endl;
    besChi2      = _chisq2sigma(besChi2);
    mixingChi2   = _chisq2sigma(mixingChi2);
    combinedChi2 = _chisq2sigma(combinedChi2);

    // Make plots of each
    std::unique_ptr<TGraph2D> besGraph =
        std::make_unique<TGraph2D>("bes", "BES", besChi2.size(), reZ.data(), imZ.data(), besChi2.data());
    std::unique_ptr<TGraph2D> mixingGraph =
        std::make_unique<TGraph2D>("mixing", "Mixing", mixingChi2.size(), reZ.data(), imZ.data(), mixingChi2.data());
    std::unique_ptr<TGraph2D> combinedGraph = std::make_unique<TGraph2D>(
        "combined", "Combination", combinedChi2.size(), reZ.data(), imZ.data(), combinedChi2.data());

    _plot(*besGraph, std::string{prefix + "bes.png"}, {reZ[minBesChi2Index], imZ[minBesChi2Index]});
    _plot(*mixingGraph, std::string{prefix + "bes_mixing.png"}, {reZ[minBesChi2Index], imZ[minBesChi2Index]});
    _plot(*combinedGraph, std::string{prefix + "bes_combined.png"}, {reZ[minBesChi2Index], imZ[minBesChi2Index]});
}

} // namespace

/*
 * Tell us what bin number to use via the CLI
 */
int main(const int argc, const char* argv[])
{
    assert(argc == 2 && "Pass 1 arg (BES phsp bin number [0-3])");

    // Generate allowed values for Z
    constexpr unsigned short N{25};
    const auto               allowedReZ = _linspace<N>(1.1);
    const auto               allowedImZ = _linspace<N>(1.1);

    // Assign nonsense values of likelihood to a special value
    constexpr double nonsense{5000.0};

    // Other parameters that we need to simulate decays etc.
    constexpr double          maxTime{0.004};
    FitterUtil::DecayParams_t params{
        .x     = CharmFitter::WORLD_AVERAGE_X,
        .y     = CharmFitter::WORLD_AVERAGE_Y,
        .r     = 0.055,
        .z_im  = std::numeric_limits<double>::quiet_NaN(),
        .z_re  = std::numeric_limits<double>::quiet_NaN(),
        .width = 2500.0,
    };

    constexpr unsigned short               numBins{4};
    std::array<const BES::Bin, numBins>    bins{BES::Bin::Bin1, BES::Bin::Bin2, BES::Bin::Bin3, BES::Bin::Bin4};
    std::array<const std::string, numBins> names{"bin0_", "bin1_", "bin2_", "bin3_"};

    const unsigned short i = std::atoi(argv[1]);
    std::cout << "Using bin " << i << std::endl;
    _makePlots(bins[i], maxTime, params, nonsense, allowedReZ, allowedImZ, names[i]);

    return 0;
}
