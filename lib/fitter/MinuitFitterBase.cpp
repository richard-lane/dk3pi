#include <iostream>

#include "D2K3PiError.h"
#include "MinuitFitterBase.h"

MinuitFitter::MinuitFitter(const FitData_t& fitData) : BaseFitter(fitData)
{
    ;
}

TMatrixD MinuitFitter::covarianceVector2CorrelationMatrix(const std::vector<double>& covarianceVector,
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

void MinuitFitter::saveFitPlot(const std::string&          plotTitle,
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

void MinuitFitter::_storeMinuitFitParams(const ROOT::Minuit2::FunctionMinimum& min)
{
    // Set our fitParams to the values obtained in the fit
    fitParams.fitParams      = min.UserParameters().Params();
    fitParams.fitParamErrors = min.UserParameters().Errors();

    // Acquire a vector representing the covariance matrix and convert it to a correlation TMatrixD
    fitParams.correlationMatrix = std::make_unique<TMatrixD>(
        covarianceVector2CorrelationMatrix(min.UserCovariance().Data(), std::vector<size_t>{}));
}