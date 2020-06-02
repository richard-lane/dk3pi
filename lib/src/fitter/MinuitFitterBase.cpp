#include <iostream>

#include "D2K3PiError.h"
#include "MinuitFitterBase.h"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"

MinuitFitterBase::MinuitFitterBase(const FitData_t& fitData) : BaseFitter(fitData)
{
    ;
}

void MinuitFitterBase::fixParameters(const std::vector<std::string>& fixParams)
{
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

    for (auto it = fixParams.begin(); it != fixParams.end(); ++it) {
        _parameters->Fix(*it);
    }
}

void MinuitFitterBase::freeParameters(const std::vector<std::string>& freeParams)
{
    size_t numFitParams  = _parameters->Parameters().size();
    size_t numFreeParams = freeParams.size();
    if (numFreeParams >= numFitParams) {
        std::cerr << "Cannot free " << numFreeParams << " parameters, only have " << numFitParams << " (";
        for (size_t i = 0; i < numFitParams; ++i) {
            std::cerr << _parameters->GetName(i) << ", ";
        }
        std::cerr << "\b\b)." << std::endl;
        throw D2K3PiException();
    }

    for (auto it = freeParams.begin(); it != freeParams.end(); ++it) {
        _parameters->Release(*it);
    }
}

void MinuitFitterBase::_setParams(const std::vector<std::string>& names,
                                  const std::vector<double>&      values,
                                  const std::vector<double>&      errors)
{
    size_t numNames  = names.size();
    size_t numValues = values.size();
    size_t numErrs   = errors.size();
    if (numNames != numValues || numErrs != numValues) {
        std::cerr << "Error setting parameters; passed " << numNames << " names, " << numValues << " values, and "
                  << numErrs << " errors" << std::endl;
        throw D2K3PiException();
    }

    _parameters = std::make_unique<ROOT::Minuit2::MnUserParameters>(values, errors);
    for (size_t i = 0; i < names.size(); i++) {
        _parameters->SetName(i, names[i]);
    }
}

void MinuitFitterBase::fit()
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

    // Create a minimiser
    ROOT::Minuit2::MnMigrad migrad(*_fitFcn, *_parameters, 2U);

    // Minimuse chi squared as defined by our _fitFcn
    min = std::make_unique<ROOT::Minuit2::FunctionMinimum>(migrad());

    // Check that our solution is "valid"
    // I think this checks that the call limit wasn't reached and that the fit converged, though it's never possible
    // to be sure with Minuit2
    if (!min->IsValid()) {
        throw BadFitException(*min);
    }

    _storeMinuitFitParams();

    // Store chi squared
    fitParams.fitStatistic = min->Fval();

    // Set our TGraph pointer to the right thing
    // The x error bars should show the bin widths
    std::vector<double> xErrors(_fitData.numBins, -1);
    for (size_t i = 0; i < _fitData.numBins; ++i) {
        xErrors[i] = 0.5 * (_fitData.binLimits[i + 1] - _fitData.binLimits[i]);
    }

    plot = std::make_unique<TGraphErrors>(
        _fitData.numBins, _fitData.binCentres.data(), _fitData.data.data(), xErrors.data(), _fitData.errors.data());
}

TMatrixD MinuitFitterBase::covarianceVector2CorrelationMatrix(const std::vector<double>& covarianceVector)
{

    if (!_parameters) {
        std::cerr << "Parameters not set" << std::endl;
        throw D2K3PiException();
    }

    // Check which parameters are fixed and create a vector of them
    std::vector<ROOT::Minuit2::MinuitParameter> minuitParams{_parameters->Trafo().Parameters()};
    std::vector<ROOT::Minuit2::MinuitParameter> fixedParams{};
    for (auto it = minuitParams.begin(); it != minuitParams.end(); ++it) {
        if (it->IsFixed()) {
            fixedParams.push_back(*it);
        }
    }

    // Check that we have the right number of elements in our covariance vector
    size_t numElements = covarianceVector.size();
    size_t numParams   = _parameters->Params().size() - fixedParams.size();

    if (!numParams) {
        std::cerr << "No free parameters found when constructing covariance matrix" << std::endl;
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
    std::vector<double> errors = fitParams.fitParamErrors;
    std::vector<size_t> paramsToRemove{};
    for (auto it = fixedParams.begin(); it != fixedParams.end(); ++it) {
        paramsToRemove.push_back(it->Number());
    }

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

    // Check that a legend has been provided
    if (!legendParams) {
        std::cerr << "Must specify legend parameters if plotting a minuit-fitted graph" << std::endl;
        throw D2K3PiException();
    }

    // Set plot title and axis labels
    bestFitFunction->SetLineColor(kRed);
    plot->SetTitle((plotTitle + ";time/ns;DCS/CF ratio").c_str());

    // Save our fit to file
    util::saveObjectsToFile<TGraphErrors>(std::vector<TObject*>{plot.get(), bestFitFunction.get()},
                                          std::vector<std::string>{"AP", "CSAME"},
                                          std::vector<std::string>{"Data", "Best fit"},
                                          path,
                                          *legendParams);
}

void MinuitFitterBase::_storeMinuitFitParams()
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
        std::make_unique<TMatrixD>(covarianceVector2CorrelationMatrix(min->UserCovariance().Data()));
}