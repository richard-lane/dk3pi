
#include <iostream>

#include "D2K3PiError.h"
#include "MinuitScannerBase.h"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnScan.h"

MinuitScannerBase::MinuitScannerBase(const FitData_t& fitData) : MinuitFitterBase(fitData)
{
    ;
}

void MinuitScannerBase::chiSqParameterScan(const size_t i, const size_t numPoints, const double low, const double high)
{
    // Check that we have sensible low and high
    if (low > high) {
        std::cerr << "Low value must be below high value for parameter scan." << std::endl;
        throw D2K3PiException();
    }

    // Check that our parameter index is in range
    if (i >= fitParams.fitParams.size()) {
        std::cerr << "Cannot scan parameter " << i << "; out of range" << std::endl;
        throw D2K3PiException();
    }

    // Check that we haven't already set _migrad to something- it will be overwritten here
    if (_migrad) {
        std::cerr << "_migrad already set- this will be overwritten by the 1d scan function." << std::endl;
        throw D2K3PiException();
    }

    // Find a vector of the parameter values we're interested in
    std::vector<double> parameterVals(numPoints, 0.0);
    double              step = (high - low) / (numPoints - 1);
    for (size_t k = 0; k < numPoints; ++k) {
        parameterVals[k] = low + k * step;
    }

    // Set parameterScan to a vector of the right length
    parameterScan = std::vector<std::pair<double, double>>(numPoints);

    // For all of these values
    for (auto k = 0; k < numPoints; ++k) {
        // Fix parameter
        _parameters->SetValue(i, parameterVals[k]);

        // Must re-set _migrad to account for the new parameter values
        _migrad = std::make_unique<ROOT::Minuit2::MnMigrad>(*_fitFcn, *_parameters);
        _migrad->Fix(i);

        // Perform a fit
        ROOT::Minuit2::FunctionMinimum min = (*_migrad)();

        // Populate chi squared
        parameterScan[k] = std::make_pair(parameterVals[k], min.Fval());

        resetMigrad();
    }
}

void MinuitScannerBase::twoDParamScan(const size_t i,
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