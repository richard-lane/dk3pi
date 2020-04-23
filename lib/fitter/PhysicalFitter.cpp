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

void PhysicalFitter::fit()
{
    // Check that parameters have been set
    if (!_parameters) {
        std::cerr << "Parameters not set" << std::endl;
        throw D2K3PiException();
    }

    // Check that we have fixed at least three of x, y, Re(Z), Im(Z) and decay width.
    // Copy the parameters we have into a vector + remove r
    std::vector<ROOT::Minuit2::MinuitParameter> x_y_rez_imz = _parameters->Trafo().Parameters();
    x_y_rez_imz.erase(x_y_rez_imz.begin() + 2);

    // Count how many of the relevant parameters are fixed
    short unsigned int numFixedParams{0};
    for (auto it = x_y_rez_imz.begin(); it != x_y_rez_imz.end(); ++it) {
        if (it->IsFixed()) {
            ++numFixedParams;
        }
    }

    // Raise if not enough are fixed
    if (numFixedParams < 3) {
        std::cerr << "Must fix at least three of {x, y, width, Im(Z), Re(Z)} for fit to be well defined ("
                  << numFixedParams << " fixed)." << std::endl;
        throw D2K3PiException();
    }

    // Use base class implementation to actually perform the fit
    MinuitFitterBase::fit();

    // Create a best-fit function
    bestFitFunction = std::make_unique<TF1>(
        "Best fit function", "[2]*[2] + [2]*([1]*[4] + [0]*[3])*[5]*x + 0.25 * ([0]*[0] + [1]*[1])*[5]*[5]*x*x");
    bestFitFunction->SetParameters(fitParams.fitParams.data());
}

void PhysicalFitter::setPhysicalFitParams(const std::vector<double>& initialParams,
                                          const std::vector<double>& initialErrors)
{
    if (initialParams.size() != 6 || initialErrors.size() != 6) {
        std::cerr << "Must provide 6 initial parameters and errors for x, y, rD, Im(Z), Re(Z) and decay width."
                  << std::endl;
        throw D2K3PiException();
    }
    _setParams(_paramNames, initialParams, initialErrors);
}
