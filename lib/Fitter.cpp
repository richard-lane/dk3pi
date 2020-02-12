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
#include "util.h"

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

void Fitter::pol2fit(const std::string& options)
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

void Fitter::expectedFunctionFit(const double minTime, const double maxTime, const std::string& options)
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
    _plot->SetTitle((plotTitle + ";time/ns;CF/DCS ratio").c_str());

    // Save our fit to file
    util::saveToFile(_plot.get(), path, "AP");
}

#endif // FITTER_CPP
