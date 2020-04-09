#include <iostream>

#include "D2K3PiError.h"
#include "MinuitPolynomialFitter.h"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"

MinuitPolynomialFitter::MinuitPolynomialFitter(const FitData_t& fitData) : MinuitFitterBase(fitData)
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
    // std::vector<double> bestFitData{_fitData.binCentres};
    // std::transform(bestFitData.begin(), bestFitData.end(), bestFitData.begin(), [&](double time) {
    //    return fitParams.fitParams[0] + fitParams.fitParams[1] * time + fitParams.fitParams[2] * time * time;
    //});
    // std::vector<double> zeros(_fitData.numPoints, 0.0); // Want errors of 0

    // bestFitPlot = std::make_unique<TGraphErrors>(
    //    _fitData.numPoints, _fitData.binCentres.data(), bestFitData.data(), zeros.data(), zeros.data());
}
