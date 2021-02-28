#include "DecaySimulator.h"
#include "fitter/ConstrainedFitter.h"
#include "physics.h"
#include "util.h"

#include "TCanvas.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH2.h"

#include "Minuit2/MnMinos.h"

/*
 * 2d scan of im(Z) and re(Z), both with and without a constraint
 */
static void test_z_scan()
{
    double                    maxTime             = 0.004;
    size_t                    numCfEvents         = 7000000;
    double                    efficiencyTimescale = 2500.0;
    FitterUtil::DecayParams_t phaseSpaceParams    = {
        .x     = CharmFitter::WORLD_AVERAGE_X,
        .y     = CharmFitter::WORLD_AVERAGE_Y,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = efficiencyTimescale,
    };

    std::random_device rd;
    std::mt19937       _gen = std::mt19937(rd());
    SimulatedDecays    MyDecays({0.0, maxTime}, phaseSpaceParams, _gen);

    // Generate data
    std::vector<double> dcsTimes{MyDecays.wsDecayTimes(Phys::numDCSDecays(numCfEvents, phaseSpaceParams, maxTime))};
    std::vector<double> cfTimes{MyDecays.rsDecayTimes(numCfEvents)};

    // Define some time bins
    std::sort(dcsTimes.begin(), dcsTimes.end());
    std::vector<double> timeBinLimits =
        util::findBinLimits(dcsTimes, 100, 0, 1.05 * maxTime); // stupid way to find bin limits

    // Create fitters
    std::array<double, 6>          initialParameterGuess{phaseSpaceParams.x,
                                                phaseSpaceParams.y,
                                                phaseSpaceParams.r,
                                                phaseSpaceParams.z_im,
                                                phaseSpaceParams.z_re,
                                                phaseSpaceParams.width};
    std::array<double, 6>          initialErrorsGuess{1, 1, 1, 1, 1, 1};
    CharmFitter::ConstrainedFitter PhysFitterConstraint(timeBinLimits, initialParameterGuess, initialErrorsGuess);

    // Fill fitter with the right data
    PhysFitterConstraint.addRSPoints(dcsTimes, std::vector(dcsTimes.size(), 1.0));
    PhysFitterConstraint.addWSPoints(cfTimes, std::vector(cfTimes.size(), 1.0));

    // No efficiency
    auto efficiency = [](const double) { return 1.0; };

    const size_t        numPoints{50};
    std::vector<double> imVals(numPoints);
    std::vector<double> reVals(numPoints);
    for (size_t i = 0; i < numPoints; ++i) {
        imVals[i] = -1.0 + (2.0 * i / numPoints);
        reVals[i] = -1.0 + (2.0 * i / numPoints);
    }

    // Perform a 2d chi squared scan on the components of Z
    PhysFitterConstraint.fixParameter("width");
    auto twoDParameterScanConstraint =
        CharmFitter::twoDParamScan(PhysFitterConstraint, efficiency, "z_im", "z_re", imVals, reVals);

    // Flatten our arrays into 1d things so we can plot them
    const size_t        numTotalPoints = numPoints * numPoints;
    std::vector<double> flatImVals(numTotalPoints);
    std::vector<double> flatReVals(numTotalPoints);
    std::vector<double> chiSquaredValsWithConstraint(numTotalPoints);
    for (size_t i = 0; i < numPoints; ++i) {
        for (size_t j = 0; j < numPoints; ++j) {
            flatImVals[i * numPoints + j]                   = imVals[i];
            flatReVals[i * numPoints + j]                   = reVals[j];
            chiSquaredValsWithConstraint[i * numPoints + j] = twoDParameterScanConstraint[i][j];
        }
    }

    // Subtract off the minimum chi squared
    double minChiSqWithConstraint =
        *std::min_element(chiSquaredValsWithConstraint.begin(), chiSquaredValsWithConstraint.end());
    std::transform(chiSquaredValsWithConstraint.begin(),
                   chiSquaredValsWithConstraint.end(),
                   chiSquaredValsWithConstraint.begin(),
                   [&](double x) { return x - minChiSqWithConstraint; });

    // Square canvas
    TCanvas* c = new TCanvas();
    c->cd();
    c->SetWindowSize(1000, 1000);

    // Draw a circle to show the allowed Z values
    double radius = 0.56; // No idea how TEllipse works but this is the factor i need to scale my ellipse by to give it
                          // a radius of 1?
    TEllipse* boundary = new TEllipse(0, 0, radius, radius);
    boundary->SetFillColorAlpha(0, 0); // Transparent

    // Take square root to find standard deviations
    std::vector<double> constrainedFitSigmaVals{chiSquaredValsWithConstraint};
    std::transform(constrainedFitSigmaVals.begin(),
                   constrainedFitSigmaVals.end(),
                   constrainedFitSigmaVals.begin(),
                   [&](double x) { return std::sqrt(x); }); // Could use std::bind

    // Draw contour plot up to 5 sigma
    size_t    numContours                = 6;
    size_t    maxSigma                   = numContours - 1;
    double    contourLevels[numContours] = {0, 1, 2, 3, 4, 5};
    TGraph2D* ConstraintGraph =
        new TGraph2D(numTotalPoints, flatImVals.data(), flatReVals.data(), constrainedFitSigmaVals.data());
    ConstraintGraph->SetMaximum(maxSigma);
    ConstraintGraph->GetHistogram()->SetContour(numContours, contourLevels);
    ConstraintGraph->SetTitle("Simulated #chi^{2} Scan of Charm Interference Parameter;Im(Z);Re(Z)");
    ConstraintGraph->GetZaxis()->SetTitle("\\sigma");

    // Point representing "true" value of Z
    TEllipse* trueZ =
        new TEllipse(radius * phaseSpaceParams.z_im, // Again have to scale by this mysterious radius factor
                     radius * phaseSpaceParams.z_re,
                     0.01,
                     0.01);

    // Save as an image
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    ConstraintGraph->Draw("CONT4Z");
    boundary->Draw();
    trueZ->Draw();
    c->SaveAs("zconstraint.png");
    delete trueZ;
    delete boundary;
    delete c;
    delete ConstraintGraph;

    // Also plot contours in the R-delta plane
    // pretty much just a copy-paste of the above code with "magphase" stuck on the end of some variable names
    std::vector<std::pair<double, double>> magAndPhase = util::reIm2magPhase(flatReVals, flatImVals);
    std::vector<double>                    magnitudes(numTotalPoints);
    std::vector<double>                    phases(numTotalPoints);
    for (size_t i = 0; i < numTotalPoints; ++i) {
        magnitudes[i] = magAndPhase[i].first;
        phases[i]     = magAndPhase[i].second;
    }

    // Square canvas
    TCanvas* c2 = new TCanvas();
    c2->cd();
    c2->SetWindowSize(1000, 1000);

    // Draw contour plot up to 5 sigma
    size_t    numContoursMagPhase                = 6;
    size_t    maxSigmaMagPhase                   = numContoursMagPhase - 1;
    double    contourLevelsMagPhase[numContours] = {0, 1, 2, 3, 4, 5};
    TGraph2D* ConstraintGraphMagPhase =
        new TGraph2D(numTotalPoints, magnitudes.data(), phases.data(), constrainedFitSigmaVals.data());
    ConstraintGraphMagPhase->SetMaximum(maxSigmaMagPhase);
    ConstraintGraphMagPhase->GetHistogram()->SetContour(numContoursMagPhase, contourLevelsMagPhase);
    ConstraintGraphMagPhase->SetTitle("Z Scan Transformed to Polar Coords;R;\\delta;\\sigmas");

    // Point representing "true" value of Z
    std::vector<std::pair<double, double>> trueZmagPhaseVals =
        util::reIm2magPhase(std::vector<double>{phaseSpaceParams.z_re}, std::vector<double>{phaseSpaceParams.z_im});
    // Need to multiply by mysterious radius factor again
    double    trueR         = radius * (2 * trueZmagPhaseVals[0].first - 1);
    double    trueDelta     = radius * (trueZmagPhaseVals[0].second / M_PI - 1);
    TEllipse* trueZMagPhase = new TEllipse(trueR, trueDelta, 0.01, 0.01);

    // Axis limits
    ConstraintGraphMagPhase->GetXaxis()->SetLimits(0, 1);

    // Save as an image
    ConstraintGraphMagPhase->Draw("CONT4Z");
    trueZMagPhase->Draw();
    c2->SaveAs("z_Rdelta.png");
}

int main()
{
    // 2d scans
    test_z_scan();
    return 0;
}
