/*
 * Small macro to test that the k3pibinning phase space gives the right value of R when using my
 * bin_generated_decays.cpp thing
 *
 * As of 14/11/19, it does
 *
 */
#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom.h"
#include "TTree.h"

#include "k3pi_binning.h"

// ---- Magic Numbers
#define NUM_EVENTS 1000

#define NUMBER_PRODUCTS 4
#define K_MASS 0.493677
#define PION_MASS 0.13957018

#define PARTICLE_MOMENTUM 0, 0, 0, 1.86962 // Decaying particle momentum in GeV
#define PARTICLE_LIFETIME 0.000425

/*
 * From a vector of TLorentzVectors, return a vector of a particle's properties
 *
 * e.g. findParticleData(i,j,myVector) will return a vector of the ith particle's j'th momentum component
 *
 * Uses the convention that 4-momentum is (px, py, pz, E)
 */
std::vector<double> findParticleData(const size_t                                    particleIndex,
                                     const size_t                                    momentumIndex,
                                     const std::vector<std::vector<TLorentzVector>> &eventsVector)
{
    size_t vectorLength = eventsVector.size();

    std::vector<double> dataVector = std::vector<double>(vectorLength);

    for (size_t i = 0; i < vectorLength; ++i) {
        dataVector[i] = eventsVector[i][particleIndex][momentumIndex];
    }

    return dataVector;
}

void example_generate_decays()
{
    // Create a phase space object to generate decays from
    TLorentzVector      dMomentum(PARTICLE_MOMENTUM);
    std::vector<double> daughterMasses = {K_MASS, PION_MASS, PION_MASS, PION_MASS};

    TGenPhaseSpace phsp;
    phsp.SetDecay(dMomentum, NUMBER_PRODUCTS, daughterMasses.data());

    // ---- Create decays using the k3pi scheme
    // Also create a random device for exponentially generated times
    std::vector<std::vector<TLorentzVector>> eventVectors = std::vector<std::vector<TLorentzVector>>(NUM_EVENTS);
    std::vector<double>                      timeVector   = std::vector<double>(NUM_EVENTS);

    std::random_device              rd;
    std::exponential_distribution<> rng(1 / PARTICLE_LIFETIME);
    std::mt19937                    rnd_gen(rd());

    for (size_t i = 0; i < NUM_EVENTS; ++i) {
        eventVectors[i] = k3pi_binning::makeUnweighted(phsp);
        timeVector[i]   = rng(rnd_gen);
    }

    // Create vectors of K and pi data
    std::vector<double> kE  = findParticleData(0, 3, eventVectors);
    std::vector<double> kPx = findParticleData(0, 0, eventVectors);
    std::vector<double> kPy = findParticleData(0, 1, eventVectors);
    std::vector<double> kPz = findParticleData(0, 2, eventVectors);

    std::vector<double> pi1E  = findParticleData(1, 3, eventVectors);
    std::vector<double> pi1Px = findParticleData(1, 0, eventVectors);
    std::vector<double> pi1Py = findParticleData(1, 1, eventVectors);
    std::vector<double> pi1Pz = findParticleData(1, 2, eventVectors);

    std::vector<double> pi2E  = findParticleData(2, 3, eventVectors);
    std::vector<double> pi2Px = findParticleData(2, 0, eventVectors);
    std::vector<double> pi2Py = findParticleData(2, 1, eventVectors);
    std::vector<double> pi2Pz = findParticleData(2, 2, eventVectors);

    std::vector<double> pi3E  = findParticleData(3, 3, eventVectors);
    std::vector<double> pi3Px = findParticleData(3, 0, eventVectors);
    std::vector<double> pi3Py = findParticleData(3, 1, eventVectors);
    std::vector<double> pi3Pz = findParticleData(3, 2, eventVectors);

    // For some reason this is what im choosing to do here
    double kEArray[NUM_EVENTS];
    double kPxArray[NUM_EVENTS];
    double kPyArray[NUM_EVENTS];
    double kPzArray[NUM_EVENTS];

    double pi1EArray[NUM_EVENTS];
    double pi1PxArray[NUM_EVENTS];
    double pi1PyArray[NUM_EVENTS];
    double pi1PzArray[NUM_EVENTS];

    double pi2EArray[NUM_EVENTS];
    double pi2PxArray[NUM_EVENTS];
    double pi2PyArray[NUM_EVENTS];
    double pi2PzArray[NUM_EVENTS];

    double pi3EArray[NUM_EVENTS];
    double pi3PxArray[NUM_EVENTS];
    double pi3PyArray[NUM_EVENTS];
    double pi3PzArray[NUM_EVENTS];

    // ---- Write them to a root file
    // This file must be initialised before the tree is created so that ROOT knows which file to write the tree to
    TFile outFile("GeneratedDecays.root", "CREATE");

    //     Create a TTree called DalitzEventList
    TTree *myTree = new TTree("DalitzEventList", "Generated D->K3pi decays");

    unsigned long long bufsize = sizeof(double) * eventVectors.size();

    for (size_t i = 0; i < kE.size(); ++i) {
        double kPxi = kPx[i];
        double kPyi = kPy[i];
        double kPzi = kPz[i];
        double kEi  = kE[i];
        myTree->Branch("_1_K~_Px", &kPxi, "kPxi/D", bufsize);
        myTree->Branch("_1_K~_Py", &kPyi, "kPyi/D", bufsize);
        myTree->Branch("_1_K~_Pz", &kPzi, "kPzi/D", bufsize);
        myTree->Branch("_1_K~_E", &kEi, "kEi/D", bufsize);

        double pi1Pxi = pi1Px[i];
        double pi1Pyi = pi1Py[i];
        double pi1Pzi = pi1Pz[i];
        double pi1Ei  = pi1E[i];
        myTree->Branch("_2_pi#_Px", &pi1Pxi, "pi1Pxi/D", bufsize);
        myTree->Branch("_2_pi#_Py", &pi1Pyi, "pi1Pyi/D", bufsize);
        myTree->Branch("_2_pi#_Pz", &pi1Pzi, "pi1Pzi/D", bufsize);
        myTree->Branch("_2_pi#_E", &pi1Ei, "pi1Ei/D", bufsize);

        double pi2Pxi = pi2Px[i];
        double pi2Pyi = pi2Py[i];
        double pi2Pzi = pi2Pz[i];
        double pi2Ei  = pi2E[i];
        myTree->Branch("_3_pi#_Px", &pi2Pxi, "pi2Pxi/D", bufsize);
        myTree->Branch("_3_pi#_Py", &pi2Pyi, "pi2Pyi/D", bufsize);
        myTree->Branch("_3_pi#_Pz", &pi2Pzi, "pi2Pzi/D", bufsize);
        myTree->Branch("_3_pi#_E", &pi2Ei, "pi2Ei/D", bufsize);

        double pi3Pxi = pi3Px[i];
        double pi3Pyi = pi3Py[i];
        double pi3Pzi = pi3Pz[i];
        double pi3Ei  = pi3E[i];
        myTree->Branch("_4_pi~_Px", &pi3Pxi, "pi3Pxi/D", bufsize);
        myTree->Branch("_4_pi~_Py", &pi3Pyi, "pi3Pyi/D", bufsize);
        myTree->Branch("_4_pi~_Pz", &pi3Pzi, "pi3Pzi/D", bufsize);
        myTree->Branch("_4_pi~_E", &pi3Ei, "pi3Ei/D", bufsize);

        double timei = timeVector[i];
        myTree->Branch("D_decayTime", &timei, "timei/D", bufsize);

        myTree->Fill();
    }

    outFile.Write();
    outFile.Close();
}
