#ifndef FITTER_CPP
#define FITTER_CPP

#include <algorithm>
#include <boost/filesystem.hpp>
#include <functional>
#include <iostream>
#include <vector>

#include "D2K3PiError.h"
#include "DecaySimulator.h"
#include "Fitter.h"
#include "MinuitFitter.h"
#include "util.h"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/VariableMetricMinimizer.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TROOT.h"

FitData::FitData() {}

FitData::FitData(const std::vector<double>& myBinCentres,
                 const std::vector<double>& myBinWidths,
                 const std::vector<double>& myData,
                 const std::vector<double>& myErrors)
{

    // Check that our bin centres are sorted
    if (!std::is_sorted(myBinCentres.begin(), myBinCentres.end())) {
        std::cerr << "Bins should be sorted" << std::endl;
        throw D2K3PiException();
    }

    // Check that all of our vectors are the same length
    size_t binCentreSize = myBinCentres.size();
    size_t binWidthSize  = myBinWidths.size();
    size_t dataLength    = myData.size();
    size_t errorsSize    = myErrors.size();
    if (binCentreSize != errorsSize || binWidthSize != errorsSize || dataLength != errorsSize) {
        std::cerr << "FitData parameters must all be the same length." << std::endl;
        throw D2K3PiException();
    }

    // Check that none of our bins overlap
    for (size_t i = 0; i < binWidthSize - 1; ++i) {
        double leftEdge  = myBinCentres[i] + 0.5 * myBinWidths[i];
        double rightEdge = myBinCentres[i + 1] - 0.5 * myBinWidths[i + 1];

        // An overlap is when the upper edge of a lower bin is higher than the lower edge of a higher bin.
        if (leftEdge > rightEdge + DBL_EPSILON) {
            std::cerr << "Bins may not overlap; bin " << i << " and bin " << i + 1 << " have values " << leftEdge
                      << " and " << rightEdge << std::endl;
            throw D2K3PiException();
        }
    }

    // Check that the data is all finite
    for (auto it = myData.begin(); it != myData.end(); ++it) {
        // Use TMath::Finite instead of std::is_finite() because for some reason including the ROOT headers makes
        // std::is_finite(-inf) sometimes return true...
        if (!TMath::Finite(*it) || *it == 0.0) {
            std::cerr << "Data must be finite: value " << *it << " encountered." << std::endl;
            throw D2K3PiException();
        }
    }

    // Now that all the checks have passed, set the class attributes
    binCentres = myBinCentres;
    data       = myData;
    errors     = myErrors;
    numPoints  = binCentreSize;

    /* Commented out as we don't want our data to have any x-uncertainty
    // Divide bin widths by 2 to get bin errors
    binErrors = myBinWidths;
    std::transform(binErrors.begin(),
                   binErrors.end(),
                   binErrors.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, 0.5));
    */
    binErrors = std::vector<double>(dataLength, 0.0);
}

Fitter::Fitter(const FitData_t& fitData)
{
    // Set our attributes to the right things
    // No need to perform consistency checks as they are performed by the FitData constructor
    _fitData = fitData;
}

void Fitter::fitUsingRootBuiltinPol2(const std::string& options)
{
    // Set our TGraph thing to the right thing
    _plot = std::make_unique<TGraphErrors>(_fitData.numPoints,
                                           _fitData.binCentres.data(),
                                           _fitData.data.data(),
                                           _fitData.binErrors.data(),
                                           _fitData.errors.data());

    // ROOT is terrible and will often fail to fit when the TGraph contains x-errors, unless the initial fit parameters
    // are reasonably close to the true values. We can get around this by perfoming an initial fit ignoring error bars
    // ("W"), then fitting again.
    TF1* graph_fit = ((TF1*)(gROOT->GetFunction("pol2")));
    _plot->Fit(graph_fit, "WQRN");

    // "S" option tells ROOT to return the result of the fit as a TFitResultPtr
    TFitResultPtr fitResult = _plot->Fit(graph_fit, (options + "S").c_str());

    // Populate the results struct
    fitParams.fitParams      = {graph_fit->GetParameter(0), graph_fit->GetParameter(1), graph_fit->GetParameter(2)};
    fitParams.fitParamErrors = {graph_fit->GetParError(0), graph_fit->GetParError(1), graph_fit->GetParError(2)};

    // Assign some memory to our correlation matrix
    fitParams.correlationMatrix = std::make_unique<TMatrixD>(fitResult->GetCorrelationMatrix());
}

void Fitter::fitUsingRootCustomFcn(const double minTime, const double maxTime, const std::string& options)
{
    // Set our TGraph pointer to the right thing
    _plot = std::make_unique<TGraphErrors>(_fitData.numPoints,
                                           _fitData.binCentres.data(),
                                           _fitData.data.data(),
                                           _fitData.binErrors.data(),
                                           _fitData.errors.data());

    // Define our custom function to fit with and initalise parameters
    // Init them all to 1 for now
    std::unique_ptr<TF1> func = std::make_unique<TF1>("func", "[0]+[1]*x+[2]*x*x", minTime, maxTime);
    func->SetParameter(0, 1.0);
    func->SetParameter(1, 1.0);
    func->SetParameter(2, 1.0);

    // ROOT is terrible and will often fail to fit when the TGraph contains x-errors, unless the initial fit parameters
    // are reasonably close to the true values. We can get around this by perfoming an initial fit ignoring error bars
    // ("W"), then fitting again.
    _plot->Fit(func.get(), "WQRN");

    // "S" option tells ROOT to return the result of the fit as a TFitResultPtr
    TFitResultPtr fitResult = _plot->Fit(func.get(), (options + "S").c_str());

    // Set fitParams to the right things
    fitParams.fitParams      = {func->GetParameter(0), func->GetParameter(1), func->GetParameter(2)};
    fitParams.fitParamErrors = {func->GetParError(0), func->GetParError(1), func->GetParError(2)};

    // Assign some memory to our correlation matrix
    fitParams.correlationMatrix = std::make_unique<TMatrixD>(fitResult->GetCorrelationMatrix());
}

TMatrixD Fitter::covarianceVector2CorrelationMatrix(const std::vector<double>& covarianceVector)
{
    // Check that we have the right number of elements in our covariance vector
    size_t numElements       = covarianceVector.size();
    size_t numFitParamErrors = fitParams.fitParamErrors.size();

    if (numFitParamErrors == 0) {
        std::cerr << "Fit has not yet been performed; fit param error vector is empty." << std::endl;
        throw D2K3PiException();
    }

    size_t expectedVectorLength = numFitParamErrors * (numFitParamErrors + 1) / 2;
    if (expectedVectorLength != numElements) {
        std::cerr << "Have " << numFitParamErrors << " fit params but " << numElements
                  << " elements in covariance vector (expected " << expectedVectorLength << ")" << std::endl;
        throw D2K3PiException();
    }

    // Create an empty TMatrixD that we will fill with the right values
    TMatrixD CorrMatrix = TMatrixD(numFitParamErrors, numFitParamErrors);

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
        double correlation      = *it / (fitParams.fitParamErrors[column] * fitParams.fitParamErrors[row]);
        CorrMatrix[column][row] = correlation;

        if (column != row) {
            CorrMatrix[row][column] = correlation;
        }
    }
    return CorrMatrix;
}

void Fitter::fitUsingMinuit2ChiSq(const std::vector<double>& initialParams, const std::vector<double>& initialErrors)
{
    // Check that we have been passed 3 initial parameters and errors
    if (initialParams.size() != 3 || initialErrors.size() != 3) {
        std::cout << "fitUsingMinuit2ChiSq requires a guess of all 3 parameters and their errors" << std::endl;
        throw D2K3PiException();
    }

    // Create an object representing our Minuit2-compatible 2nd order polynomial
    BasePolynomialFcn FitFcn(_fitData.data, _fitData.binCentres, _fitData.errors);

    // Create a minimiser and minimise our chi squared
    ROOT::Minuit2::VariableMetricMinimizer Minimizer;
    ROOT::Minuit2::FunctionMinimum         min = Minimizer.Minimize(FitFcn, initialParams, initialErrors);

    // Check that our solution is "valid"
    // I think this checks that the call limit wasn't reached and that the fit converged, though it's never possible to
    // be sure with Minuit2
    if (!min.IsValid()) {
        std::cerr << "Minuit fit invalid" << std::endl;
        // It would be nice to print out the Minimiser status here but std::cout << doesn't work (the docs are lying)
        throw D2K3PiException();
    }

    // Store our fit parameters and correlation matrix as class attributes
    _storeMinuitFitParams(min);

    // Set our TGraph pointer to the right thing
    _plot = std::make_unique<TGraphErrors>(_fitData.numPoints,
                                           _fitData.binCentres.data(),
                                           _fitData.data.data(),
                                           _fitData.binErrors.data(),
                                           _fitData.errors.data());

    // Create also a best-fit dataset from our parameters and data, plotting this on the same
    std::vector<double> bestFitData{_fitData.binCentres};
    std::transform(bestFitData.begin(), bestFitData.end(), bestFitData.begin(), [&](double time) {
        return fitParams.fitParams[0] + fitParams.fitParams[1] * time + fitParams.fitParams[2] * time * time;
    });
    std::vector<double> zeros(_fitData.numPoints, 0.0); // Want errors of 0

    _bestFitPlot = std::make_unique<TGraphErrors>(
        _fitData.numPoints, _fitData.binCentres.data(), bestFitData.data(), zeros.data(), zeros.data());
}

void Fitter::_storeMinuitFitParams(const ROOT::Minuit2::FunctionMinimum& min)
{

    // Set our fitParams to the values obtained in the fit
    // I couldn't find an obvious way to do this so... here we are
    fitParams.fitParams = std::vector<double>{min.UserParameters().Parameter(0).Value(),
                                              min.UserParameters().Parameter(1).Value(),
                                              min.UserParameters().Parameter(2).Value()};

    fitParams.fitParamErrors = std::vector<double>{min.UserParameters().Parameter(0).Error(),
                                                   min.UserParameters().Parameter(1).Error(),
                                                   min.UserParameters().Parameter(2).Error()};

    // Acquire a vector representing the covariance matrix and convert it to a correlation TMatrixD
    fitParams.correlationMatrix =
        std::make_unique<TMatrixD>(covarianceVector2CorrelationMatrix(min.UserCovariance().Data()));
}

void Fitter::saveFitPlot(const std::string& plotTitle, const std::string& path)
{
    // Check that a fit has been made
    if (fitParams.fitParams.empty() || _plot == nullptr) {
        std::cerr << "Run the fitter before plotting" << std::endl;
        throw D2K3PiException();
    }

    // Check that path doesn't already exist
    if (boost::filesystem::exists(path)) {
        std::cerr << path << " already exists" << std::endl;
        throw D2K3PiException();
    }

    // Set the plot title; for some reason this is also how to set axis labels?
    _plot->SetTitle((plotTitle + ";time/ns;DCS/CF ratio").c_str());

    // Save our fit to file
    // If root's builtin was used, _bestFitPlot will not have been assigned
    if (_bestFitPlot != nullptr) {
        _bestFitPlot->SetLineColor(kRed);
        util::saveObjectsToFile<TGraphErrors>(
            std::vector<TObject*>{_plot.get(), _bestFitPlot.get()}, std::vector<std::string>{"AP", "CSAME"}, path);

    } else {
        util::saveObjectToFile(_plot.get(), path, "AP");
    }
}

#endif // FITTER_CPP
