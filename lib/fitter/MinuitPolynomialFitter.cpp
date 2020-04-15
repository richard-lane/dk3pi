#include <iostream>

#include "D2K3PiError.h"
#include "MinuitPolynomialFitter.h"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "TF1.h"

MinuitPolynomialFitter::MinuitPolynomialFitter(const FitData_t& fitData) : MinuitScannerBase(fitData)
{
    _fitFcn = std::make_unique<PolynomialChiSqFcn>(_fitData.data, _fitData.binCentres, _fitData.errors);
}

void MinuitPolynomialFitter::fit(const std::vector<size_t>& fixParams)
{
    if (!_parameters) {
        std::cerr << "parameters not set; call this->setPolynomialParams before fitting." << std::endl;
        throw D2K3PiException();
    }

    MinuitFitterBase::fit(fixParams);

    // Set our TGraph pointer to the right thing
    plot = std::make_unique<TGraphErrors>(_fitData.numPoints,
                                          _fitData.binCentres.data(),
                                          _fitData.data.data(),
                                          _fitData.binErrors.data(),
                                          _fitData.errors.data());

    bestFitFunction = std::make_unique<TF1>("Best fit function", "[0] +[1]*x+[2]*x*x");
    bestFitFunction->SetParameters(fitParams.fitParams.data());
}

void MinuitPolynomialFitter::setPolynomialParams(const std::vector<double>& initialParams,
                                                 const std::vector<double>& initialErrors)
{
    _setParams(std::vector<std::string>{"a", "b", "c"}, initialParams, initialErrors);
}
