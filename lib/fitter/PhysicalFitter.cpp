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

    // Check that we have fixed at least one of x, y, Re(Z) or Im(Z)- otherwise our fit is poorly defined
    // Create a vector of Minuit's internal representations of x, y, Re(Z) and Im(Z)
    std::vector<ROOT::Minuit2::MinuitParameter> x_y_rez_imz = _parameters->Trafo().Parameters();
    x_y_rez_imz.erase(x_y_rez_imz.begin() + 2);
    x_y_rez_imz.erase(x_y_rez_imz.begin() + 5);

    // Check that at least one of them is fixed
    if (std::all_of(x_y_rez_imz.begin(), x_y_rez_imz.end(), [](ROOT::Minuit2::MinuitParameter& param) {
            return !param.IsFixed();
        })) {
        std::cerr << "Must fix one of x, y, or a component of Z for fit to be well defined" << std::endl;
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
    _setParams(_paramNames, initialParams, initialErrors);
}