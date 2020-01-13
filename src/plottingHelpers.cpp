/*
 * Helper functions for plotting to check consistency when binning
 *
 */

#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLorentzVector.h"

/*
 * From a vector of TLorentzVectors and the desired index (0,1,2,3), find a C-style array of data
 *
 * Allocates memory to the array which must be freed by the caller.
 *
 * e.g. vector2Array(myVector, 0) for x-momentum
 */
double *vector2Array(const std::vector<TLorentzVector> &particleVector, const size_t index)
{
    size_t  length   = particleVector.size();
    double *outArray = new double[length];

    for (size_t i = 0; i < length; ++i) {
        outArray[i] = particleVector[i][index];
    }

    return outArray;
}

/*
 * CoM energy squared of the ab system
 * Takes two vectors of TLorentzVectors as particle data
 *
 * Assumes each entry of the vector is of the form (Px, Py, Pz, E)
 *
 */
const std::vector<double> s(const std::vector<TLorentzVector> &particleA, const std::vector<TLorentzVector> &particleB)
{
    size_t              length = particleA.size();
    std::vector<double> sValues(particleA.size());

    for (size_t i = 0; i < length; ++i) {

        sValues[i] = std::pow(particleA[i][3], 2) - std::pow(particleA[i][0], 2) - std::pow(particleA[i][1], 2) -
                     std::pow(particleA[i][2], 2) + std::pow(particleB[i][3], 2) - std::pow(particleB[i][0], 2) -
                     std::pow(particleB[i][1], 2) - std::pow(particleB[i][2], 2) +
                     2 * particleA[i][3] * particleB[i][3] - 2 * particleA[i][0] * particleB[i][0] -
                     2 * particleA[i][1] * particleB[i][1] - 2 * particleA[i][2] * particleB[i][2];
    }

    return sValues;
}

/*
 * Plot an n-bin histogram from an array
 */
void plot_hist(const std::string &title,
               double *           myData,
               size_t             length,
               const float        xmin,
               const float        xmax,
               const int          nBins)
{

    const char *titleStr = title.c_str();
    auto        kCanvas  = new TCanvas(titleStr, titleStr, 600, 600);
    TH1D *      hist     = new TH1D(titleStr, titleStr, nBins, xmin, xmax);
    hist->FillN(length, myData, 0);
    hist->Draw();
}

/*
 * Plot the K energies and make a plot of s01 vs s02 to check consistency with the ROOT TBrowser
 *
 */
void plot_things(const std::vector<TLorentzVector> &kVectors,
                 const std::vector<TLorentzVector> &pi1Vectors,
                 const std::vector<TLorentzVector> &pi2Vectors)
{
    size_t length = kVectors.size();

    // Plot K energies
    plot_hist("K energies", vector2Array(kVectors, 3), length, 0.45, 1, 100);

    // Plot CoM energies on a new canvas
    auto                      comCanvas = new TCanvas("CoM Energies", "CoM Energies", 600, 600);
    const std::vector<double> s01       = s(kVectors, pi1Vectors);
    const std::vector<double> s02       = s(kVectors, pi2Vectors);
    TGraph *                  myGraph   = new TGraph(length, s01.data(), s02.data());
    myGraph->Draw("AP");
}
