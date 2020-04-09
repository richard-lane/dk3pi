#include <iostream>

#include "D2K3PiError.h"
#include "MinuitPolynomialFitter.h"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"

MinuitPolynomialFitter::MinuitPolynomialFitter(const FitData_t& fitData) : MinuitFitterBase(fitData)
{
    ;
}

void MinuitPolynomialFitter::fit(const std::vector<double>&                    initialParams,
                                 const std::vector<double>&                    initialErrors,
                                 const FitAlgorithm_t&                         FitMethod,
                                 const std::vector<std::pair<size_t, double>>& fixParams)
{
    // Check that we have been passed 3 initial parameters and errors
    if (initialParams.size() != 3 || initialErrors.size() != 3) {
        std::cout << "fitUsingMinuit2ChiSq requires a guess of all 3 parameters and their errors" << std::endl;
        throw D2K3PiException();
    }

    // Create an object representing our Minuit2-compatible 2nd order polynomial
    _fitFcn = std::make_unique<PolynomialChiSqFcn>(_fitData.data, _fitData.binCentres, _fitData.errors);

    MinuitFitterBase::fit(initialParams, initialErrors, FitMethod, fixParams);

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
