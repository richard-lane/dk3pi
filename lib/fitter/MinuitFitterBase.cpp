#include <iostream>

#include "D2K3PiError.h"
#include "MinuitFitterBase.h"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"

MinuitFitterBase::MinuitFitterBase(const FitData_t& fitData) : BaseFitter(fitData)
{
    ;
}

void MinuitFitterBase::setParams(const std::vector<std::string>& names,
                                 const std::vector<double>&      values,
                                 const std::vector<double>&      errors)
{
    _parameters = std::make_unique<ROOT::Minuit2::MnUserParameters>(values, errors);
    for (size_t i = 0; i < names.size(); i++) {
        _parameters->SetName(i, names[i]);
    }
}

void MinuitFitterBase::fit(const std::vector<size_t>& fixParams)
{
    if (!_fitFcn) {
        std::cerr << "Fit fcn not yet set - cannot perform fit. Classes inheriting from MinuitFitterBase must set "
                     "the _fitFcn attribute."
                  << std::endl;
        throw D2K3PiException();
    }

    if (!_parameters) {
        std::cerr << "Parameters not yet set- have you called setParams()?" << std::endl;
        throw D2K3PiException();
    }

    // Check that we aren't fixing more params than we have
    size_t numFitParams = _parameters->Parameters().size();
    size_t numFixParams = fixParams.size();
    if (numFixParams >= numFitParams) {
        std::cerr << "Cannot fix " << numFixParams << " parameters, only have " << numFitParams << " (";
        for (size_t i = 0; i < numFitParams; ++i) {
            std::cerr << _parameters->GetName(i) << ", ";
        }
        std::cerr << "\b\b)." << std::endl;
        throw D2K3PiException();
    }

    // Create a minimiser
    ROOT::Minuit2::MnMigrad migrad(*_fitFcn, *_parameters);

    // If we find index i in our list of parameters to fix, fix it.
    for (auto it = fixParams.begin(); it != fixParams.end(); ++it) {
        migrad.Fix(*it);
    }

    // Minimuse chi squared as defined by our _fitFcn
    min = std::make_unique<ROOT::Minuit2::FunctionMinimum>(migrad());

    // Check that our solution is "valid"
    // I think this checks that the call limit wasn't reached and that the fit converged, though it's never possible
    // to be sure with Minuit2
    if (!min->IsValid()) {
        std::cerr << "Minuit fit invalid" << std::endl;
        std::cerr << *min << std::endl;
        throw D2K3PiException();
    }

    _storeMinuitFitParams(fixParams);

    // Store chi squared
    statistic = std::make_unique<double>(min->Fval());

    // Set our TGraph pointer to the right thing
    plot = std::make_unique<TGraphErrors>(_fitData.numPoints,
                                          _fitData.binCentres.data(),
                                          _fitData.data.data(),
                                          _fitData.binErrors.data(),
                                          _fitData.errors.data());
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
    if (!legendParams) {
        std::cerr << "Must specify legend parameters if plotting a minuit-fitted graph" << std::endl;
        throw D2K3PiException();
    }
    // bestFitPlot->SetLineColor(kRed);
    // util::saveObjectsToFile<TGraphErrors>(std::vector<TObject*>{plot.get(), bestFitPlot.get()},
    //                                      std::vector<std::string>{"AP", "CSAME"},
    //                                      std::vector<std::string>{"Data", "Best fit"},
    //                                      path,
    //                                      *legendParams);
}

void MinuitFitterBase::_storeMinuitFitParams(const std::vector<size_t>& fixParamIndices)
{
    if (!min) {
        std::cerr << "Cannot store params: fit not yet run" << std::endl;
        throw D2K3PiException();
    }

    // Set our fitParams to the values obtained in the fit
    fitParams.fitParams      = min->UserParameters().Params();
    fitParams.fitParamErrors = min->UserParameters().Errors();

    // Acquire a vector representing the covariance matrix and convert it to a correlation TMatrixD
    fitParams.correlationMatrix =
        std::make_unique<TMatrixD>(covarianceVector2CorrelationMatrix(min->UserCovariance().Data(), fixParamIndices));
}