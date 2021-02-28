/*
 * Plot a likelihood scan using the CLEO constraint only
 */
#include "ConstrainedFitter.h"
#include "cleo_interface.h"
#include "util.h"

#include <TGraph2D.h>
#include <TH2D.h>

int main(void)
{
    // Create arrays covering the range of Re(Z) and Im(Z) values
    constexpr unsigned short  N{100}; // Actually we make N+1 points
    std::array<double, N + 1> allowedReZ{};
    std::array<double, N + 1> allowedImZ{};

    for (size_t i = 0; i < N + 1; ++i) {
        double zVal   = -1.0 + (i * 2.0) / N;
        allowedReZ[i] = zVal;
        allowedImZ[i] = zVal;
    }

    // Easier to store our points in three arrays
    std::vector<double> reZ{};
    std::vector<double> imZ{};
    std::vector<double> likelihood{};

    // Assign nonsense values to a likelihood of 0.0
    constexpr double nonsense = 0.0;

    // For each pair, find the CLEO likelihood
    for (const auto r : allowedReZ) {
        for (const auto i : allowedImZ) {
            reZ.push_back(r);
            imZ.push_back(i);

            // Construct some phsp params to use
            FitterUtil::DecayParams_t phaseSpaceParams{
                .x     = CharmFitter::WORLD_AVERAGE_X,
                .y     = CharmFitter::WORLD_AVERAGE_Y,
                .r     = 0.055,
                .z_im  = i,
                .z_re  = r,
                .width = 2500.0,
            };

            double l = -2 * CLEO::cleoLikelihood(CLEO::Bin::Bin1, phaseSpaceParams);
            l        = std::isnan(l) ? nonsense : l;

            likelihood.push_back(l);
        }
    }

    // Find min and max likelihoods
    double min = 1e6;
    double max = 0.0;
    for (const auto l : likelihood) {
        if (l == nonsense)
            continue;

        if (l < min) {
            min = l;
        }
        if (l > max) {
            max = l;
        }
    }

    // Remove the min likelihood from each value
    std::transform(likelihood.begin(), likelihood.end(), likelihood.begin(), [&min](double x) { return x - min; });

    // Take the sqrt to transform into std dev
    std::transform(likelihood.begin(), likelihood.end(), likelihood.begin(), [&](double x) {
        if (x < 0.0)
            return x;
        return std::sqrt(x);
    }); // Could use std::bind

    std::unique_ptr<TH1D> hist = std::make_unique<TH1D>("", "", 100, -5.0, 260.0);
    hist->FillN(likelihood.size(), likelihood.data(), nullptr);
    util::saveObjectToFile(hist.get(), "hist.png");

    // Create a contour plot of the CLEO likelihood
    constexpr short int             numContours{7};
    std::array<double, numContours> contours{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    std::unique_ptr<TGraph2D> graph =
        std::make_unique<TGraph2D>("contour", "", likelihood.size(), reZ.data(), imZ.data(), likelihood.data());

    // Hack
    graph->SetMinimum(0.0);
    graph->SetMaximum(std::sqrt(max));

    graph->GetHistogram()->SetContour(numContours, contours.data());
    graph->SetTitle("CLEO constraint only (Phsp Bin 0)");
    graph->GetZaxis()->SetTitle("\\sigma");
    graph->GetXaxis()->SetTitle("Re(Z)");
    graph->GetYaxis()->SetTitle("Im(Z)");

    std::unique_ptr<TCanvas> c = std::make_unique<TCanvas>();
    c->cd();
    c->SetWindowSize(1000, 1000);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.15);
    graph->Draw("CONT4Z");
    c->SaveAs("cleo.png");

    return 0;
}
