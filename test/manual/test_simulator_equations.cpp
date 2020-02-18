/*
 * Run a simulation using a SimulatedDecays object and draw a graph of the two curves.
 */
#include <utility>

#include <TH1D.h>

#include "DecaySimulator.h"
#include "util.h"

/*
 * Set the decay parameters to sensible values.
 */
void setParams(DecayParams_t &DecayParams)
{
    // Global Params
    DecayParams.x     = 0.004;
    DecayParams.y     = 0.007;
    DecayParams.width = 2500.0; // Width is in nanoseconds

    // Per-bin params
    DecayParams.r    = 1.14;
    DecayParams.z_re = 0.7;
    DecayParams.z_im = -0.3;
}

/*
 * Create an object representing the RS and WS decays and plot their time dependence using sensible values for the
 * particle and phase-space parameters.
 */
void simulateDecays()
{
    // Create a struct of sensible parameter values to use.
    DecayParams_t MyParams;
    setParams(MyParams);

    // Allow times up to maxTime nanoseconds and rates up to 1.5 (arbitary units; we only care about ratios)
    double                    maxTime      = 0.002;
    std::pair<double, double> allowedTimes = std::make_pair(0, maxTime);
    std::pair<double, double> allowedRates = std::make_pair(0, 1.3);

    // Create our Decay simulator object
    size_t          numDecays = 10000;
    SimulatedDecays MyDecays  = SimulatedDecays(allowedTimes, allowedRates, MyParams);

    MyDecays.plotRates();

    // Keep generating points, checking if they are in either distribution until both the RS and WS vectors are
    // populated.
    MyDecays.findCfDecayTimes(numDecays);
    MyDecays.findDcsDecayTimes(numDecays);

    // Do some stuff with histograms
    size_t numBins = 20;
    TH1D * RSHist  = new TH1D("Test accept-reject, RS", "", numBins, 0, maxTime);
    TH1D * WSHist  = new TH1D("Test accept-reject, WS", "", numBins, 0, maxTime);

    for (size_t i = 0; i < numDecays; ++i) {
        RSHist->Fill(MyDecays.RSDecayTimes[i]);
        WSHist->Fill(MyDecays.WSDecayTimes[i]);
    }

    util::saveToFile(RSHist, "RSHist.pdf");
    util::saveToFile(WSHist, "WSHist.pdf");
    delete RSHist;
    delete WSHist;
}

int main()
{
    simulateDecays();

    return 0;
}
