/*
 * Plot a chi2 scan using the CLEO + BES-III constraints
 */
#include "ConstrainedFitter.h"
#include "bes_interface.h"
#include "charm_interface.h"
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
    std::vector<double> chi2{};

    // Assign nonsense values to the current minimum chi2
    double nonsense = 1e6;

    // Just so i only have to change the code in 1 place
    const short                           binNumber{4}; // 1, 2, 3, or 4
    std::map<const short, const BES::Bin> bins{
        {1, BES::Bin::Bin1}, {2, BES::Bin::Bin2}, {3, BES::Bin::Bin3}, {4, BES::Bin::Bin4}};

    // For each pair, find the bes chi2
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

            double c = charm_threshhold::chi2(bins[binNumber], phaseSpaceParams);
            if (c < nonsense)
                nonsense = c;

            // If this likelihood value is NaN, set it to the current minimum
            if (std::isnan(c))
                c = nonsense;

            chi2.push_back(c);
        }
    }

    // Remove the min chi2 from each value
    const auto min = *std::min_element(chi2.begin(), chi2.end());
    std::cout << min << std::endl;
    std::transform(chi2.begin(), chi2.end(), chi2.begin(), [&min](double x) { return x - min; });
    std::cout << min << std::endl;

    std::unique_ptr<TH1D> hist = std::make_unique<TH1D>("", "", 100, -5.0, 260.0);
    hist->FillN(chi2.size(), chi2.data(), nullptr);
    hist->SetTitle("Charm Threshhold chi2;chi2;Count");
    util::saveObjectToFile(hist.get(), (std::string{"charm_hist_"} + std::to_string(binNumber) + ".png").c_str());

    // Create a contour plot of the BES chi2
    constexpr short int             numContours{7};
    std::array<double, numContours> contours{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    std::unique_ptr<TGraph2D> graph =
        std::make_unique<TGraph2D>("contour", "", chi2.size(), reZ.data(), imZ.data(), chi2.data());

    // Hack
    const auto max = *std::max_element(chi2.begin(), chi2.end());
    graph->SetMinimum(0.0);
    graph->SetMaximum(max);

    graph->GetHistogram()->SetContour(numContours, contours.data());
    graph->SetTitle((std::string{"BES + CLEO constraint (Phsp Bin "} + std::to_string(binNumber) + ")").c_str());
    graph->GetZaxis()->SetTitle("\\sigma");
    graph->GetXaxis()->SetTitle("Re(Z)");
    graph->GetYaxis()->SetTitle("Im(Z)");

    std::unique_ptr<TCanvas> c = std::make_unique<TCanvas>();
    c->cd();
    c->SetWindowSize(1000, 1000);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.15);
    graph->Draw("CONT4Z");
    c->SaveAs((std::string{"charm_"} + std::to_string(binNumber) + ".png").c_str());

    return 0;
}
