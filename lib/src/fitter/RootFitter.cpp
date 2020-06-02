#include "RootFitter.h"
#include "D2K3PiError.h"
#include "util.h"

#include "boost/filesystem.hpp"

#include "TF1.h"
#include "TFitResult.h"

RootFitter::RootFitter(const FitData_t& fitData) : BaseFitter(fitData)
{
    ;
}

void RootFitter::fit(const double minTime, const double maxTime, const std::string& options)
{
    // Set our TGraph pointer to the right thing
    // The x error bars should show the bin widths
    std::vector<double> xErrors(_fitData.numBins, -1);
    for (size_t i = 0; i < _fitData.numBins; ++i) {
        xErrors[i] = 0.5 * (_fitData.binLimits[i + 1] - _fitData.binLimits[i]);
    }
    plot = std::make_unique<TGraphErrors>(
        _fitData.numBins, _fitData.binCentres.data(), _fitData.data.data(), xErrors.data(), _fitData.errors.data());

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
    fitParams.fitStatistic = plot->Chisquare(func.get());
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