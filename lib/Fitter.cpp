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
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnUserParameters.h"
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

BaseFitter::BaseFitter(const FitData_t& fitData)
{
    // Set our attributes to the right things
    // No need to perform consistency checks as they are performed by the FitData constructor
    _fitData = fitData;
}

RootFitter::RootFitter(const FitData_t& fitData) : BaseFitter(fitData)
{
    ;
}

void RootFitter::fit(const double minTime, const double maxTime, const std::string& options)
{
    // Set our TGraph pointer to the right thing
    plot = std::make_unique<TGraphErrors>(_fitData.numPoints,
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
    plot->Fit(func.get(), "WQRN");

    // "S" option tells ROOT to return the result of the fit as a TFitResultPtr
    TFitResultPtr fitResult = plot->Fit(func.get(), (options + "S").c_str());

    // Set fitParams to the right things
    fitParams.fitParams      = {func->GetParameter(0), func->GetParameter(1), func->GetParameter(2)};
    fitParams.fitParamErrors = {func->GetParError(0), func->GetParError(1), func->GetParError(2)};

    // Assign some memory to our correlation matrix
    fitParams.correlationMatrix = std::make_unique<TMatrixD>(fitResult->GetCorrelationMatrix());

    // Set chi squared to the value stored in func
    statistic = std::make_unique<double>(plot->Chisquare(func.get()));
}

void RootFitter::saveFitPlot(const std::string& plotTitle, const std::string& path)
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

    util::saveObjectToFile(plot.get(), path, "AP");
}

MinuitPolynomialFitter::MinuitPolynomialFitter(const FitData_t& fitData) : BaseFitter(fitData)
{
    ;
}

void MinuitPolynomialFitter::fit(const std::vector<double>& initialParams,
                                 const std::vector<double>& initialErrors,
                                 const FitAlgorithm_t&      FitMethod)
{
    // Check that we have been passed 3 initial parameters and errors
    if (initialParams.size() != 3 || initialErrors.size() != 3) {
        std::cout << "fitUsingMinuit2ChiSq requires a guess of all 3 parameters and their errors" << std::endl;
        throw D2K3PiException();
    }

    // Create an object representing our Minuit2-compatible 2nd order polynomial
    _fitFcn = std::make_unique<PolynomialChiSqFcn>(_fitData.data, _fitData.binCentres, _fitData.errors);

    // Store our parameters
    ROOT::Minuit2::MnUserParameters parameters;
    parameters.Add("a", initialParams[0], initialErrors[0]);
    parameters.Add("b", initialParams[1], initialErrors[1]);
    parameters.Add("c", initialParams[2], initialErrors[2]);

    // Create a minimiser and minimise our chi squared
    ROOT::Minuit2::MnMigrad migrad(*_fitFcn, parameters);
    min = std::make_unique<ROOT::Minuit2::FunctionMinimum>(migrad());

    // Check that our solution is "valid"
    // I think this checks that the call limit wasn't reached and that the fit converged, though it's never possible to
    // be sure with Minuit2
    if (!min->IsValid()) {
        std::cerr << "Minuit fit invalid" << std::endl;
        std::cerr << *min << std::endl;
        throw D2K3PiException();
    }

    // Store our fit parameters and correlation matrix as class attributes
    _storeMinuitFitParams(*min);

    // Store chi squared
    statistic = std::make_unique<double>(min->Fval());

    // Set our TGraph pointer to the right thing
    plot = std::make_unique<TGraphErrors>(_fitData.numPoints,
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

    bestFitPlot = std::make_unique<TGraphErrors>(
        _fitData.numPoints, _fitData.binCentres.data(), bestFitData.data(), zeros.data(), zeros.data());
}

TMatrixD MinuitPolynomialFitter::covarianceVector2CorrelationMatrix(const std::vector<double>& covarianceVector)
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

void MinuitPolynomialFitter::saveFitPlot(const std::string&          plotTitle,
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

void MinuitPolynomialFitter::_storeMinuitFitParams(const ROOT::Minuit2::FunctionMinimum& min)
{
    // Set our fitParams to the values obtained in the fit
    fitParams.fitParams      = min.UserParameters().Params();
    fitParams.fitParamErrors = min.UserParameters().Errors();

    // Acquire a vector representing the covariance matrix and convert it to a correlation TMatrixD
    fitParams.correlationMatrix =
        std::make_unique<TMatrixD>(covarianceVector2CorrelationMatrix(min.UserCovariance().Data()));
}

MinuitPolyScan::MinuitPolyScan(const FitData_t& fitData) : MinuitPolynomialFitter(fitData)
{
    ;
}

void MinuitPolyScan::chiSqParameterScan(const size_t i, const size_t numPoints, const double low, const double high)
{
    // Check that a fit has been performed
    if (!_fitFcn) {
        std::cerr << "Must run Minuit fitter before performing parameter scan." << std::endl;
        throw D2K3PiException();
    }

    // Check that we have sensible low and high
    if (low > high) {
        std::cerr << "Low value must be below high value for parameter scan." << std::endl;
        throw D2K3PiException();
    }

    // Check that we don't have too many points
    if (numPoints > 100) {
        std::cerr << "Cannot scan a parameter at more than 100 points (sorry)." << std::endl;
        throw D2K3PiException();
    }

    // Check that our parameter index is in range
    if (i >= fitParams.fitParams.size()) {
        std::cerr << "Cannot scan parameter " << i << "; out of range" << std::endl;
        throw D2K3PiException();
    }

    // Create MnScan object that will be used to perform our scan and run the scan
    ROOT::Minuit2::MnScan Scanner = ROOT::Minuit2::MnScan(*_fitFcn, fitParams.fitParams, fitParams.fitParamErrors);
    parameterScan                 = Scanner.Scan(i, numPoints + 1, low, high);
}

void MinuitPolyScan::twoDParamScan(const size_t i,
                                   const size_t j,
                                   const size_t iPoints,
                                   const size_t jPoints,
                                   const double iLow,
                                   const double iHigh,
                                   const double jLow,
                                   const double jHigh)
{
    if (i >= fitParams.fitParams.size() || j >= fitParams.fitParams.size()) {
        std::cerr << "Cannot scan params " << i << ", " << j << "; only have " << fitParams.fitParams.size()
                  << " params." << std::endl;
        throw D2K3PiException();
    }

    if (i == j) {
        std::cerr << "Cannot perfom 2d scan on a parameter against itself" << std::endl;
        throw D2K3PiException();
    }

    if (!_fitFcn) {
        std::cerr << "Run fit before running 2d parameter scan" << std::endl;
        throw D2K3PiException();
    }

    // Create vectors to store the i, j values we're interested in, and a vector to store our result in
    std::vector<double> iVals(iPoints);
    std::vector<double> jVals(jPoints);
    twoDParameterScan = std::vector<std::vector<double>>(iPoints * jPoints);

    // Fill vectors of i and j values
    // Technically we could do all this in the big loop where we do the minimisation but i think this is clearer to
    // think about
    double iStep = (iHigh - iLow) / (iPoints - 1);
    double jStep = (jHigh - jLow) / (jPoints - 1);
    for (size_t k = 0; k < iPoints; ++k) {
        iVals[k] = iLow + k * iStep;
    }
    for (size_t k = 0; k < jPoints; ++k) {
        jVals[k] = jLow + k * jStep;
    }

    // Loop over our desired i and j values, perform a fit with i and j fixed + fill twoDParameterScan
    for (size_t iPoint = 0; iPoint < iPoints; ++iPoint) {
        for (size_t jPoint = 0; jPoint < jPoints; ++jPoint) {
            double iVal = iVals[iPoint];
            double jVal = jVals[jPoint];

            // Minimise wrt the params
            ROOT::Minuit2::MnUserParameters parameters;
            for (size_t param = 0; param < fitParams.fitParams.size(); ++param) {
                parameters.Add(std::to_string(param), fitParams.fitParams[param], fitParams.fitParamErrors[param]);
            }

            parameters.SetValue(i, iVal);
            parameters.SetValue(j, jVal);

            // Create a minimiser
            ROOT::Minuit2::MnMigrad migrad(*_fitFcn, parameters);

            // Hack: if we have 6 params, assume one is width + should be fixed
            // TODO get rid of
            if (fitParams.fitParams.size() == 6) {
                migrad.Fix(5);
            }

            // Fix i and j; now they are not allowed to vary from the value they get set to
            migrad.Fix(i);
            migrad.Fix(j);

            ROOT::Minuit2::FunctionMinimum min = migrad();
            if (!min.IsValid()) {
                std::cerr << "Minuit fit invalid" << std::endl;
                std::cout << min << std::endl;
                throw D2K3PiException();
            }

            twoDParameterScan[jPoint + jPoints * iPoint] = std::vector<double>{iVal, jVal, min.Fval()};
        }
    }
}

PhysicalFitter::PhysicalFitter(const FitData_t& fitData) : BaseFitter(fitData)
{
    ;
}

void PhysicalFitter::fit(const std::vector<double>& initialParams,
                         const std::vector<double>& initialErrors,
                         const FitAlgorithm_t&      FitMethod,
                         const std::vector<size_t>  fixParams)
{
    std::cout << "This isnt implemented properly yet" << std::endl;
    throw D2K3PiException();

    // Check that we have been passed 6 initial parameters and errors
    if (initialParams.size() != 6 || initialErrors.size() != 6) {
        std::cout << "fit requires a guess of 6 parameters and their errors" << std::endl;
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

    // Create a minimiser and minimise our chi squared
    ROOT::Minuit2::MnMigrad migrad(*_fitFcn, parameters);
    migrad.Fix(5); // Fix width
    ROOT::Minuit2::FunctionMinimum min = migrad();

    // Check that our solution is "valid"
    // I think this checks that the call limit wasn't reached and that the fit converged, though it's never possible to
    // be sure with Minuit2
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
    // Hack: remove width from fit params, store matrix and add it again
    double widthErr = fitParams.fitParamErrors[5];
    fitParams.fitParamErrors.pop_back();
    fitParams.correlationMatrix =
        std::make_unique<TMatrixD>(covarianceVector2CorrelationMatrix(min.UserCovariance().Data()));
    fitParams.fitParamErrors.push_back(widthErr);

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

void PhysicalFitter::saveFitPlot(const std::string&          plotTitle,
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
    bestFitPlot->SetLineColor(kRed);
    util::saveObjectsToFile<TGraphErrors>(std::vector<TObject*>{plot.get(), bestFitPlot.get()},
                                          std::vector<std::string>{"AP", "CSAME"},
                                          std::vector<std::string>{"Data", "Best fit"},
                                          path,
                                          *legendParams);
}

TMatrixD PhysicalFitter::covarianceVector2CorrelationMatrix(const std::vector<double>& covarianceVector)
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

void PhysicalFitter::_storeMinuitFitParams(const ROOT::Minuit2::FunctionMinimum& min)
{
    // Set our fitParams to the values obtained in the fit
    fitParams.fitParams      = min.UserParameters().Params();
    fitParams.fitParamErrors = min.UserParameters().Errors();

    // Acquire a vector representing the covariance matrix and convert it to a correlation TMatrixD
    fitParams.correlationMatrix =
        std::make_unique<TMatrixD>(covarianceVector2CorrelationMatrix(min.UserCovariance().Data()));
}

ParamScanner::ParamScanner(const FitData_t& fitData) : PhysicalFitter(fitData)
{
    ;
}

void ParamScanner::chiSqParameterScan(const size_t i, const size_t numPoints, const double low, const double high)
{
    // Check that a fit has been performed
    if (!_fitFcn) {
        std::cerr << "Must run Minuit fitter before performing parameter scan." << std::endl;
        throw D2K3PiException();
    }

    // Check that we have sensible low and high
    if (low > high) {
        std::cerr << "Low value must be below high value for parameter scan." << std::endl;
        throw D2K3PiException();
    }

    // Check that we don't have too many points
    if (numPoints > 100) {
        std::cerr << "Cannot scan a parameter at more than 100 points (sorry)." << std::endl;
        throw D2K3PiException();
    }

    // Check that our parameter index is in range
    if (i >= fitParams.fitParams.size()) {
        std::cerr << "Cannot scan parameter " << i << "; out of range" << std::endl;
        throw D2K3PiException();
    }

    // Create MnScan object that will be used to perform our scan and run the scan
    ROOT::Minuit2::MnScan Scanner = ROOT::Minuit2::MnScan(*_fitFcn, fitParams.fitParams, fitParams.fitParamErrors);
    parameterScan                 = Scanner.Scan(i, numPoints + 1, low, high);
}

void ParamScanner::twoDParamScan(const size_t i,
                                 const size_t j,
                                 const size_t iPoints,
                                 const size_t jPoints,
                                 const double iLow,
                                 const double iHigh,
                                 const double jLow,
                                 const double jHigh)
{
    if (i >= fitParams.fitParams.size() || j >= fitParams.fitParams.size()) {
        std::cerr << "Cannot scan params " << i << ", " << j << "; only have " << fitParams.fitParams.size()
                  << " params." << std::endl;
        throw D2K3PiException();
    }

    if (i == j) {
        std::cerr << "Cannot perfom 2d scan on a parameter against itself" << std::endl;
        throw D2K3PiException();
    }

    if (!_fitFcn) {
        std::cerr << "Run fit before running 2d parameter scan" << std::endl;
        throw D2K3PiException();
    }

    // Create vectors to store the i, j values we're interested in, and a vector to store our result in
    std::vector<double> iVals(iPoints);
    std::vector<double> jVals(jPoints);
    twoDParameterScan = std::vector<std::vector<double>>(iPoints * jPoints);

    // Fill vectors of i and j values
    // Technically we could do all this in the big loop where we do the minimisation but i think this is clearer to
    // think about
    double iStep = (iHigh - iLow) / (iPoints - 1);
    double jStep = (jHigh - jLow) / (jPoints - 1);
    for (size_t k = 0; k < iPoints; ++k) {
        iVals[k] = iLow + k * iStep;
    }
    for (size_t k = 0; k < jPoints; ++k) {
        jVals[k] = jLow + k * jStep;
    }

    // Loop over our desired i and j values, perform a fit with i and j fixed + fill twoDParameterScan
    for (size_t iPoint = 0; iPoint < iPoints; ++iPoint) {
        for (size_t jPoint = 0; jPoint < jPoints; ++jPoint) {
            double iVal = iVals[iPoint];
            double jVal = jVals[jPoint];

            // Minimise wrt the params
            ROOT::Minuit2::MnUserParameters parameters;
            for (size_t param = 0; param < fitParams.fitParams.size(); ++param) {
                parameters.Add(std::to_string(param), fitParams.fitParams[param], fitParams.fitParamErrors[param]);
            }

            parameters.SetValue(i, iVal);
            parameters.SetValue(j, jVal);

            // Create a minimiser
            ROOT::Minuit2::MnMigrad migrad(*_fitFcn, parameters);

            // Hack: if we have 6 params, assume one is width + should be fixed
            // TODO get rid of
            if (fitParams.fitParams.size() == 6) {
                migrad.Fix(5);
            }

            // Fix i and j; now they are not allowed to vary from the value they get set to
            migrad.Fix(i);
            migrad.Fix(j);

            ROOT::Minuit2::FunctionMinimum min = migrad();
            if (!min.IsValid()) {
                std::cerr << "Minuit fit invalid" << std::endl;
                std::cout << min << std::endl;
                throw D2K3PiException();
            }

            twoDParameterScan[jPoint + jPoints * iPoint] = std::vector<double>{iVal, jVal, min.Fval()};
        }
    }
}

#endif // FITTER_CPP
