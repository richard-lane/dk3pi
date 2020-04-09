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
#include "MinuitFitter.h"
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

MinuitPolynomialFitter::MinuitPolynomialFitter(const FitData_t& fitData) : MinuitFitter(fitData)
{
    ;
}

void MinuitPolynomialFitter::fit(const std::vector<double>& initialParams,
                                 const std::vector<double>& initialErrors,
                                 const FitAlgorithm_t&      FitMethod)
{
    // Check that we have been passed 3 initial parameters and errors
    if (initialParams.size() != 3 || initialErrors.size() != 3) {
        std::cout << "fitUsingMinuit2ChiSq requires a guess of all 3 parameters and their errors" << std::endl;
        throw D2K3PiException();
    }

    // Create an object representing our Minuit2-compatible 2nd order polynomial
    _fitFcn = std::make_unique<PolynomialChiSqFcn>(_fitData.data, _fitData.binCentres, _fitData.errors);

    // Store our parameters
    ROOT::Minuit2::MnUserParameters parameters;
    parameters.Add("a", initialParams[0], initialErrors[0]);
    parameters.Add("b", initialParams[1], initialErrors[1]);
    parameters.Add("c", initialParams[2], initialErrors[2]);

    // Create a minimiser and minimise our chi squared
    ROOT::Minuit2::MnMigrad migrad(*_fitFcn, parameters);
    min = std::make_unique<ROOT::Minuit2::FunctionMinimum>(migrad());

    // Check that our solution is "valid"
    // I think this checks that the call limit wasn't reached and that the fit converged, though it's never possible to
    // be sure with Minuit2
    if (!min->IsValid()) {
        std::cerr << "Minuit fit invalid" << std::endl;
        std::cerr << *min << std::endl;
        throw D2K3PiException();
    }

    // Store our fit parameters and correlation matrix as class attributes
    _storeMinuitFitParams(*min);

    // Store chi squared
    statistic = std::make_unique<double>(min->Fval());

    // Set our TGraph pointer to the right thing
    plot = std::make_unique<TGraphErrors>(_fitData.numPoints,
                                          _fitData.binCentres.data(),
                                          _fitData.data.data(),
                                          _fitData.binErrors.data(),
                                          _fitData.errors.data());

    // Create also a best-fit dataset from our parameters and data, plotting this on the same
    std::vector<double> bestFitData{_fitData.binCentres};
    std::transform(bestFitData.begin(), bestFitData.end(), bestFitData.begin(), [&](double time) {
        return fitParams.fitParams[0] + fitParams.fitParams[1] * time + fitParams.fitParams[2] * time * time;
    });
    std::vector<double> zeros(_fitData.numPoints, 0.0); // Want errors of 0

    bestFitPlot = std::make_unique<TGraphErrors>(
        _fitData.numPoints, _fitData.binCentres.data(), bestFitData.data(), zeros.data(), zeros.data());
}

MinuitPolyScan::MinuitPolyScan(const FitData_t& fitData) : MinuitPolynomialFitter(fitData)
{
    ;
}

void MinuitPolyScan::chiSqParameterScan(const size_t i, const size_t numPoints, const double low, const double high)
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

void MinuitPolyScan::twoDParamScan(const size_t i,
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

PhysicalFitter::PhysicalFitter(const FitData_t& fitData) : MinuitFitter(fitData)
{
    ;
}

void PhysicalFitter::fit(const std::vector<double>&                    initialParams,
                         const std::vector<double>&                    initialErrors,
                         const FitAlgorithm_t&                         FitMethod,
                         const std::vector<std::pair<size_t, double>>& fixParams)
{
    // Check that we have been passed 6 initial parameters and errors
    if (initialParams.size() != 6 || initialErrors.size() != 6) {
        std::cout << "fit requires a guess of 6 parameters and their errors" << std::endl;
        throw D2K3PiException();
    }

    if (fixParams.empty()) {
        std::cerr << "Must fix one of x, y, Re(Z) or Im(Z) when fitting with rD, x, y etc." << std::endl;
        throw D2K3PiException();
    }

    size_t numFixParams = fixParams.size();
    if (numFixParams > 5) {
        std::cerr << "cannot fix more than 5 parameters" << std::endl;
        throw D2K3PiException();
    }

    // Check that we have fixed at least one of x, y, Re(Z) or Im(Z)- otherwise our fit is poorly defined
    std::vector<size_t> fixParamIndices(fixParams.size());
    std::vector<double> fixParamValues(fixParams.size());
    for (size_t i = 0; i < fixParams.size(); ++i) {
        fixParamIndices[i] = fixParams[i].first;
        fixParamValues[i]  = fixParams[i].second;
    }
    if (std::find(fixParamIndices.begin(), fixParamIndices.end(), 0) == fixParamIndices.end() && // x
        std::find(fixParamIndices.begin(), fixParamIndices.end(), 1) == fixParamIndices.end() && // y
        std::find(fixParamIndices.begin(), fixParamIndices.end(), 3) == fixParamIndices.end() && // z_im
        std::find(fixParamIndices.begin(), fixParamIndices.end(), 4) == fixParamIndices.end()    // z_re
    ) {
        std::cerr << "Must fix one of x, y, or a component of Z for fit to be well defined" << std::endl;
        throw D2K3PiException();
    }

    // Create an object representing our Minuit2-compatible 2nd order polynomial
    _fitFcn = std::make_unique<DetailedPolynomialChiSqFcn>(_fitData.data, _fitData.binCentres, _fitData.errors);

    // Store our parameters
    ROOT::Minuit2::MnUserParameters parameters;
    parameters.Add("x", initialParams[0], initialErrors[0]);
    parameters.Add("y", initialParams[1], initialErrors[1]);
    parameters.Add("r", initialParams[2], initialErrors[2]);
    parameters.Add("z_im", initialParams[3], initialErrors[3]);
    parameters.Add("z_re", initialParams[4], initialErrors[4]);
    parameters.Add("width", initialParams[5], initialErrors[5]);

    // Create a minimiser and fix any parameters to the right values
    ROOT::Minuit2::MnMigrad migrad(*_fitFcn, parameters);
    for (auto it = fixParams.begin(); it != fixParams.end(); ++it) {
        parameters.SetValue(it->first, it->second);
        migrad.Fix(it->first);
    }

    // Minimuse chi squared as defined by our _fitFcn
    ROOT::Minuit2::FunctionMinimum min = migrad();

    // Check that our solution is "valid"
    // I think this checks that the call limit wasn't reached and that the fit converged, though it's never possible
    // to be sure with Minuit2
    if (!min.IsValid()) {
        std::cerr << "Minuit fit invalid" << std::endl;
        std::cerr << min << std::endl;
        throw D2K3PiException();
    }

    // Store our fit parameters and correlation matrix as class attributes
    // Set our fitParams to the values obtained in the fit
    fitParams.fitParams      = min.UserParameters().Params();
    fitParams.fitParamErrors = min.UserParameters().Errors();

    // Acquire a vector representing the covariance matrix and convert it to a correlation TMatrixD
    fitParams.correlationMatrix =
        std::make_unique<TMatrixD>(covarianceVector2CorrelationMatrix(min.UserCovariance().Data(), fixParamIndices));

    // Store chi squared
    statistic = std::make_unique<double>(min.Fval());

    // Set our TGraph pointer to the right thing
    plot = std::make_unique<TGraphErrors>(_fitData.numPoints,
                                          _fitData.binCentres.data(),
                                          _fitData.data.data(),
                                          _fitData.binErrors.data(),
                                          _fitData.errors.data());

    // Create also a best-fit dataset from our parameters and data, plotting this on the same
    DecayParams_t bestFitParams = DecayParameters{.x     = fitParams.fitParams[0],
                                                  .y     = fitParams.fitParams[1],
                                                  .r     = fitParams.fitParams[2],
                                                  .z_im  = fitParams.fitParams[3],
                                                  .z_re  = fitParams.fitParams[4],
                                                  .width = fitParams.fitParams[5]};

    // Should use std::transform
    std::vector<double> bestFitData(_fitData.binCentres.size());
    std::vector<double> zeros(_fitData.numPoints, 0.0); // Want errors of 0
    for (size_t i = 0; i < bestFitData.size(); ++i) {
        bestFitData[i] = fitPolynomial(bestFitParams, _fitData.binCentres[i]);
    }

    bestFitPlot = std::make_unique<TGraphErrors>(
        _fitData.numPoints, _fitData.binCentres.data(), bestFitData.data(), zeros.data(), zeros.data());
}

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
