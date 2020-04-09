#include <iostream>

#include "D2K3PiError.h"
#include "MinuitFitterBase.h"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"

MinuitFitterBase::MinuitFitterBase(const FitData_t& fitData) : BaseFitter(fitData)
{
    ;
}

void MinuitFitterBase::fit(const std::vector<double>&                    initialParams,
                           const std::vector<double>&                    initialErrors,
                           const FitAlgorithm_t&                         FitMethod,
                           const std::vector<std::pair<size_t, double>>& fixParams)
{
    size_t numFixParams = fixParams.size();
    if (numFixParams > 5) {
        std::cerr << "cannot fix more than 5 parameters" << std::endl;
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
    std::vector<size_t> fixParamIndices(fixParams.size());
    for (size_t i = 0; i < fixParams.size(); ++i) {
        fixParamIndices[i] = fixParams[i].first;
    }
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

TMatrixD MinuitFitterBase::covarianceVector2CorrelationMatrix(const std::vector<double>& covarianceVector,
                                                              const std::vector<size_t>& fixedParams)
{
    // Check that we have the right number of elements in our covariance vector
    size_t numElements = covarianceVector.size();
    size_t numParams   = fitParams.fitParamErrors.size() - fixedParams.size();

    if (numParams == 0) {
        std::cerr << "Fit has not yet been performed; fit param error vector is empty (or the fit has been run with 0 "
                     "free parameters)."
                  << std::endl;
        throw D2K3PiException();
    }

    size_t expectedVectorLength = numParams * (numParams + 1) / 2;
    if (expectedVectorLength != numElements) {
        std::cerr << "Have " << numParams << " fit params but " << numElements
                  << " elements in covariance vector (expected " << expectedVectorLength << ")" << std::endl;
        throw D2K3PiException();
    }

    // Create an empty TMatrixD that we will fill with the right values
    TMatrixD CorrMatrix = TMatrixD(numParams, numParams);

    // Find which error values are relevant; i.e. those which correspond to free parameters
    // Do this by copying the vector + removing the values in order, highest-first
    // Highest-first to avoid issues with looping over a vector + deleting elements from it concurrently.
    std::vector<double> errors         = fitParams.fitParamErrors;
    std::vector<size_t> paramsToRemove = fixedParams;
    std::sort(paramsToRemove.rbegin(), paramsToRemove.rend());
    for (auto it = paramsToRemove.begin(); it != paramsToRemove.end(); ++it) {
        errors.erase(errors.begin() + *it);
    }

    // We need to divide each element in our vector of covariances with the standard deviation of two parameters
    // We will need to find the position in the matrix of each element in our covariance vector, so we know which
    // errors to divide by
    size_t column = -1; // We just want a number such that when we add 1 we get 0; unsigned int overflow is safe!
    size_t row    = 0;
    for (auto it = covarianceVector.begin(); it != covarianceVector.end(); ++it) {
        column++;
        if (column > row) {
            row++;
            column = 0;
        }
        // Now that we know which variances to divide by, let's do it
        double correlation      = *it / (errors[column] * errors[row]);
        CorrMatrix[column][row] = correlation;

        if (column != row) {
            CorrMatrix[row][column] = correlation;
        }
    }
    return CorrMatrix;
}

void MinuitFitterBase::saveFitPlot(const std::string&          plotTitle,
                                   const std::string&          path,
                                   const util::LegendParams_t* legendParams)
{
    // Check that a fit has been made
    if (fitParams.fitParams.empty() || plot == nullptr) {
        std::cerr << "Run the fitter before plotting" << std::endl;
        throw D2K3PiException();
    }

    // Check that path doesn't already exist
    if (boost::filesystem::exists(path)) {
        std::cerr << path << " already exists" << std::endl;
        throw D2K3PiException();
    }

    // Set the plot title; for some reason this is also how to set axis labels?
    plot->SetTitle((plotTitle + ";time/ns;DCS/CF ratio").c_str());

    // Save our fit to file
    // If root's builtin was used, _bestFitPlot will not have been assigned
    if (bestFitPlot != nullptr) {
        if (!legendParams) {
            std::cerr << "Must specify legend parameters if plotting a minuit-fitted graph" << std::endl;
            throw D2K3PiException();
        }
        bestFitPlot->SetLineColor(kRed);
        util::saveObjectsToFile<TGraphErrors>(std::vector<TObject*>{plot.get(), bestFitPlot.get()},
                                              std::vector<std::string>{"AP", "CSAME"},
                                              std::vector<std::string>{"Data", "Best fit"},
                                              path,
                                              *legendParams);

    } else {
        util::saveObjectToFile(plot.get(), path, "AP");
    }
}

void MinuitFitterBase::_storeMinuitFitParams(const ROOT::Minuit2::FunctionMinimum& min)
{
    // Set our fitParams to the values obtained in the fit
    fitParams.fitParams      = min.UserParameters().Params();
    fitParams.fitParamErrors = min.UserParameters().Errors();

    // Acquire a vector representing the covariance matrix and convert it to a correlation TMatrixD
    fitParams.correlationMatrix = std::make_unique<TMatrixD>(
        covarianceVector2CorrelationMatrix(min.UserCovariance().Data(), std::vector<size_t>{}));
}