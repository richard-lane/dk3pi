#ifndef FITTER_CPP
#define FITTER_CPP

#include <algorithm>
#include <boost/filesystem.hpp>
#include <functional>
#include <iostream>
#include <vector>

#include "D2K3PiError.h"
#include "DecaySimulator.h"
#include "Fitter.h"
#include "FitterUtils.h"
#include "MinuitFitter.h"
#include "MinuitFitterBase.h"
#include "MinuitPolynomialFitter.h"
#include "PhysicalFitter.h"
#include "util.h"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/VariableMetricMinimizer.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TROOT.h"

ParamScanner::ParamScanner(const FitData_t& fitData) : PhysicalFitter(fitData)
{
    ;
}

void ParamScanner::chiSqParameterScan(const size_t i, const size_t numPoints, const double low, const double high)
{
    // Check that a fit has been performed
    if (!_fitFcn) {
        std::cerr << "Must run Minuit fitter before performing parameter scan." << std::endl;
        throw D2K3PiException();
    }

    // Check that we have sensible low and high
    if (low > high) {
        std::cerr << "Low value must be below high value for parameter scan." << std::endl;
        throw D2K3PiException();
    }

    // Check that we don't have too many points
    if (numPoints > 100) {
        std::cerr << "Cannot scan a parameter at more than 100 points (sorry)." << std::endl;
        throw D2K3PiException();
    }

    // Check that our parameter index is in range
    if (i >= fitParams.fitParams.size()) {
        std::cerr << "Cannot scan parameter " << i << "; out of range" << std::endl;
        throw D2K3PiException();
    }

    // Create MnScan object that will be used to perform our scan and run the scan
    ROOT::Minuit2::MnScan Scanner = ROOT::Minuit2::MnScan(*_fitFcn, fitParams.fitParams, fitParams.fitParamErrors);
    parameterScan                 = Scanner.Scan(i, numPoints + 1, low, high);
}

void ParamScanner::twoDParamScan(const size_t i,
                                 const size_t j,
                                 const size_t iPoints,
                                 const size_t jPoints,
                                 const double iLow,
                                 const double iHigh,
                                 const double jLow,
                                 const double jHigh)
{
    if (i >= fitParams.fitParams.size() || j >= fitParams.fitParams.size()) {
        std::cerr << "Cannot scan params " << i << ", " << j << "; only have " << fitParams.fitParams.size()
                  << " params." << std::endl;
        throw D2K3PiException();
    }

    if (i == j) {
        std::cerr << "Cannot perfom 2d scan on a parameter against itself" << std::endl;
        throw D2K3PiException();
    }

    if (!_fitFcn) {
        std::cerr << "Run fit before running 2d parameter scan" << std::endl;
        throw D2K3PiException();
    }

    // Create vectors to store the i, j values we're interested in, and a vector to store our result in
    std::vector<double> iVals(iPoints);
    std::vector<double> jVals(jPoints);
    twoDParameterScan = std::vector<std::vector<double>>(iPoints * jPoints);

    // Fill vectors of i and j values
    // Technically we could do all this in the big loop where we do the minimisation but i think this is clearer to
    // think about
    double iStep = (iHigh - iLow) / (iPoints - 1);
    double jStep = (jHigh - jLow) / (jPoints - 1);
    for (size_t k = 0; k < iPoints; ++k) {
        iVals[k] = iLow + k * iStep;
    }
    for (size_t k = 0; k < jPoints; ++k) {
        jVals[k] = jLow + k * jStep;
    }

    // Loop over our desired i and j values, perform a fit with i and j fixed + fill twoDParameterScan
    for (size_t iPoint = 0; iPoint < iPoints; ++iPoint) {
        for (size_t jPoint = 0; jPoint < jPoints; ++jPoint) {
            double iVal = iVals[iPoint];
            double jVal = jVals[jPoint];

            // Minimise wrt the params
            ROOT::Minuit2::MnUserParameters parameters;
            for (size_t param = 0; param < fitParams.fitParams.size(); ++param) {
                parameters.Add(std::to_string(param), fitParams.fitParams[param], fitParams.fitParamErrors[param]);
            }

            parameters.SetValue(i, iVal);
            parameters.SetValue(j, jVal);

            // Create a minimiser
            ROOT::Minuit2::MnMigrad migrad(*_fitFcn, parameters);

            // Hack: if we have 6 params, assume one is width + should be fixed
            // TODO get rid of
            if (fitParams.fitParams.size() == 6) {
                migrad.Fix(5);
            }

            // Fix i and j; now they are not allowed to vary from the value they get set to
            migrad.Fix(i);
            migrad.Fix(j);

            ROOT::Minuit2::FunctionMinimum min = migrad();
            if (!min.IsValid()) {
                std::cerr << "Minuit fit invalid" << std::endl;
                std::cout << min << std::endl;
                throw D2K3PiException();
            }

            twoDParameterScan[jPoint + jPoints * iPoint] = std::vector<double>{iVal, jVal, min.Fval()};
        }
    }
}

#endif // FITTER_CPP
