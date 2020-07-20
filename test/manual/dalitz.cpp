#include <random>

#include <TFile.h>
#include <TGenPhaseSpace.h>
#include <TH2D.h>

#include "efficiencyUtil.h"
#include "flatPhsp.h"
#include "util.h"

/*
 * Generate 3body events
 */
int main()
{
    double                parentMass = 3;
    std::array<double, 3> dauaghterMasses{0.5, 0.5, 0.5};

    // Set up phase space to be X -> B B B
    TLorentzVector stationaryParent(0, 0, 0, parentMass);
    TGenPhaseSpace phaseSpace;
    phaseSpace.SetDecay(stationaryParent, dauaghterMasses.size(), dauaghterMasses.data());

    // Create some random number generators and stuff
    std::random_device                     rd;
    std::mt19937                           gen(rd());
    std::uniform_real_distribution<double> uniformDistribution(0, phaseSpace.GetWtMax());

    // Generate a load of events
    size_t numEvents = 50000;

    // Represent an event as a vector of kinematic params; multiple events are a vector of these vectors
    std::vector<std::vector<kinematicParams_t>> events(numEvents);
    std::vector<double>                         weights(numEvents);
    for (size_t i = 0; i < numEvents; ++i) {
        events[i] = randomEvent(phaseSpace, &gen, uniformDistribution, &weights[i]);
    }

    // Calculate invariant masses for these events
    std::vector<double> m12(numEvents);
    std::vector<double> m23(numEvents);
    std::vector<double> m31(numEvents);
    for (size_t i = 0; i < numEvents; ++i) {
        m12[i] = invariantMass(std::vector<kinematicParams_t>{events[i][0], events[i][1]});
        m23[i] = invariantMass(std::vector<kinematicParams_t>{events[i][1], events[i][2]});
    }

    // Create and fill histograms for these events (Dalitz plots)
    size_t numBins    = 100;
    double min        = 0.5;
    double max        = 3;
    TH2D*  m12m23Hist = new TH2D("m12", "m12 vs m23", numBins, min, max, numBins, min, max);
    m12m23Hist->FillN(numEvents, m12.data(), m23.data(), nullptr);

    //  Save Dalitz plots to file
    util::saveObjectToFile(m12m23Hist, "m12m23.png", "COLZ");

    delete m12m23Hist;
    return 0;
}
