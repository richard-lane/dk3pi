#include <iostream>

#include "D2K3PiError.h"
#include "PhysicalFitter.h"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "TF1.h"

PhysicalFitter::PhysicalFitter(const FitData_t& fitData) : MinuitScannerBase(fitData)
{
    _fitFcn = std::make_unique<DetailedPolynomialChiSqFcn>(_fitData.data, _fitData.binCentres, _fitData.errors);
}

void PhysicalFitter::fit(const std::vector<size_t>& fixParams)
{
    // Check that parameters have been set
    if (!_parameters) {
        std::cerr << "Parameters not set" << std::endl;
        throw D2K3PiException();
    }

    // Check that we have fixed at least one of x, y, Re(Z) or Im(Z)- otherwise our fit is poorly defined
    if (fixParams.empty() || (std::find(fixParams.begin(), fixParams.end(), 0) == fixParams.end() && // x
                              std::find(fixParams.begin(), fixParams.end(), 1) == fixParams.end() && // y
                              std::find(fixParams.begin(), fixParams.end(), 3) == fixParams.end() && // z_im
                              std::find(fixParams.begin(), fixParams.end(), 4) == fixParams.end())   // z_re
    ) {
        std::cerr << "Must fix one of x, y, or a component of Z for fit to be well defined" << std::endl;
        throw D2K3PiException();
    }

    // Use base class implementation to actually perform the fit
    // For now always fix width
    std::vector<size_t> allFixParams{fixParams};
    if (std::find(fixParams.begin(), fixParams.end(), 5) == fixParams.end()) {
        allFixParams.push_back(5);
    }
    MinuitFitterBase::fit(allFixParams);

    // Create a best-fit function
    bestFitFunction = std::make_unique<TF1>(
        "Best fit function", "[2]*[2] + [2]*([1]*[4] + [0]*[3])*[5]*x + 0.25 * ([0]*[0] + [1]*[1])*[5]*[5]*x*x");
    bestFitFunction->SetParameters(fitParams.fitParams.data());
}

void PhysicalFitter::setPhysicalFitParams(const std::vector<double>& initialParams,
                                          const std::vector<double>& initialErrors)
{
    _setParams(_paramNames, initialParams, initialErrors);
}