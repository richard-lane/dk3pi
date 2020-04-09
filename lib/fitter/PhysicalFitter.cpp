#include <iostream>

#include "D2K3PiError.h"
#include "PhysicalFitter.h"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"

PhysicalFitter::PhysicalFitter(const FitData_t& fitData) : MinuitFitterBase(fitData)
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

    // Check that we have fixed at least one of x, y, Re(Z) or Im(Z)- otherwise our fit is poorly defined
    std::vector<size_t> fixParamIndices(fixParams.size());
    std::vector<double> fixParamValues(fixParams.size());
    for (size_t i = 0; i < fixParams.size(); ++i) {
        fixParamIndices[i] = fixParams[i].first;
        fixParamValues[i]  = fixParams[i].second;
    }
    if (fixParams.empty() ||
        (std::find(fixParamIndices.begin(), fixParamIndices.end(), 0) == fixParamIndices.end() && // x
         std::find(fixParamIndices.begin(), fixParamIndices.end(), 1) == fixParamIndices.end() && // y
         std::find(fixParamIndices.begin(), fixParamIndices.end(), 3) == fixParamIndices.end() && // z_im
         std::find(fixParamIndices.begin(), fixParamIndices.end(), 4) == fixParamIndices.end())   // z_re
    ) {
        std::cerr << "Must fix one of x, y, or a component of Z for fit to be well defined" << std::endl;
        throw D2K3PiException();
    }

    // Create an object representing our Minuit2-compatible 2nd order polynomial
    // TODO move to constructor or something
    _fitFcn = std::make_unique<DetailedPolynomialChiSqFcn>(_fitData.data, _fitData.binCentres, _fitData.errors);

    // Use base class implementation to actually perform the fit
    MinuitFitterBase::fit(initialParams, initialErrors, FitMethod, fixParams);

    // Create also a best-fit dataset from our parameters and data, plotting this on the same
    // DecayParams_t bestFitParams = DecayParameters{.x     = fitParams.fitParams[0],
    //                                              .y     = fitParams.fitParams[1],
    //                                              .r     = fitParams.fitParams[2],
    //                                              .z_im  = fitParams.fitParams[3],
    //                                              .z_re  = fitParams.fitParams[4],
    //                                              .width = fitParams.fitParams[5]};

    //// Should use std::transform
    // std::vector<double> bestFitData(_fitData.binCentres.size());
    // std::vector<double> zeros(_fitData.numPoints, 0.0); // Want errors of 0
    // for (size_t i = 0; i < bestFitData.size(); ++i) {
    //    bestFitData[i] = fitPolynomial(bestFitParams, _fitData.binCentres[i]);
    //}

    // bestFitPlot = std::make_unique<TGraphErrors>(
    //    _fitData.numPoints, _fitData.binCentres.data(), bestFitData.data(), zeros.data(), zeros.data());
}