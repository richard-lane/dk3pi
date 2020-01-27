#ifndef FITTER_CPP
#define FITTER_CPP

#include <algorithm>
#include <boost/filesystem.hpp>
#include <iostream>
#include <vector>

#include "../include/D2K3PiError.h"
#include "../include/Fitter.h"
#include "../include/util.h"

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
    // This might cause issues when bins share edges so I might have to think of something Cleverer
    for (size_t i = 0; i < binWidthSize - 1; ++i) {
        double leftEdge  = myBinCentres[i] + 0.5 * myBinWidths[i];
        double rightEdge = myBinCentres[i + 1] - 0.5 * myBinWidths[i + 1];
        if (leftEdge > rightEdge) {
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
    binWidths  = myBinWidths;
    data       = myData;
    errors     = myErrors;
    numPoints  = binCentreSize;
}

Fitter::Fitter(const FitData_t& fitData)
{
    // Set our attributes to the right things
    // No need to perform consistency checks as they are performed by the FitData constructor
    _fitData = fitData;
}

void Fitter::pol2fit(void)
{
    // Set our TGraph thing to the right thing
    _plot = TGraphErrors(_fitData.numPoints,
                         _fitData.binCentres.data(),
                         _fitData.data.data(),
                         _fitData.binWidths.data(),
                         _fitData.errors.data());

    // ROOT is terrible and will often fail to fit when the TGraph contains x-errors, unless the initial fit parameters
    // are reasonably close to the true values. We can get around this by perfoming an initial fit ignoring error bars
    // ("W"), then fitting again.
    TF1* graph_fit = ((TF1*)(gROOT->GetFunction("pol2")));
    _plot.Fit(graph_fit, "WQRN");

    // "S" option tells ROOT to return the result of the fit as a TFitResultPtr
    TFitResultPtr fitResult = _plot.Fit(graph_fit, "S");

    // Populate the results struct
    fitParams.fitParams      = {graph_fit->GetParameter(0), graph_fit->GetParameter(1), graph_fit->GetParameter(2)};
    fitParams.fitParamErrors = {graph_fit->GetParError(0), graph_fit->GetParError(1), graph_fit->GetParError(2)};

    // TODO fix this, it doesn't work at the moment
    fitParams.correlationMatrix = fitResult->GetCorrelationMatrix();
}

void Fitter::saveFitPlot(const std::string& plotTitle, const std::string& path)
{
    // Check that a fit has been made
    if (fitParams.fitParams.empty()) {
        std::cerr << "Run the fitter before plotting" << std::endl;
        throw D2K3PiException();
    }

    // Check that path doesn't already exist
    if (boost::filesystem::exists(path)) {
        std::cerr << path << " already exists" << std::endl;
        throw D2K3PiException();
    }

    // Set the plot title; for some reason this is also how to set axis labels?
    _plot.SetTitle((plotTitle + ";time/ns;CF/DCS ratio").c_str());

    // Save our fit to file
    util::saveToFile(&_plot, path, "AP");
}

#endif // FITTER_CPP
