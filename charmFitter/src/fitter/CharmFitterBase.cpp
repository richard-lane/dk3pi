#include <TF1.h>

#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include "CharmFitterBase.h"
#include "FitError.h"
#include "util.h"

namespace CharmFitter
{
CharmBaseFcn::CharmBaseFcn(const std::vector<double>& data,
                           const std::vector<double>& times,
                           const std::vector<double>& errors,
                           const IntegralOptions_t&   integralOptions)
    : theMeasurements(data), thePositions(times), theMVariances(errors), theErrorDef(1.),
      _integralOptions(integralOptions)
{
    if (_integralOptions.binLimits.size() != thePositions.size() + 1) {
        std::cerr << _integralOptions.binLimits.size() << " elements to bin limits but " << thePositions.size()
                  << " datapoints" << std::endl;
        throw D2K3PiException{};
    }
}

CharmFitterBase::CharmFitterBase(const std::vector<double>& binLimits)
    : _numBins{binLimits.size() - 1}, _binEdges(binLimits)
{
    // Check our bin limits are sorted
    if (!std::is_sorted(binLimits.begin(), binLimits.end())) {
        std::cerr << "Bins should be sorted" << std::endl;
        throw D2K3PiException();
    }

    // Initialise vectors of counts to the right lengths
    _rsCounts = std::vector<double>(_numBins, 0.0);
    _wsCounts = std::vector<double>(_numBins, 0.0);

    // Populate bin widths + bin centres
    _binWidths  = std::vector<double>(_numBins);
    _binCentres = std::vector<double>(_numBins);

    for (size_t i = 0; i < _numBins; ++i) {
        _binWidths[i]  = _binEdges[i + 1] - _binEdges[i];
        _binCentres[i] = (_binEdges[i + 1] + _binEdges[i]) / 2;
    }
}

void CharmFitterBase::addRSPoint(const double time, const double weight)
{
    _rsCounts[util::findBinIndex(_binEdges, time)] += weight;
}

void CharmFitterBase::addWSPoint(const double time, const double weight)
{
    _wsCounts[util::findBinIndex(_binEdges, time)] += weight;
}

void CharmFitterBase::fixParameter(const std::string& paramName)
{
    _parameters->Fix(paramName);
}

void CharmFitterBase::freeParameter(const std::string& paramName)
{
    _parameters->Release(paramName);
}

void CharmFitterBase::_setParams(const std::vector<std::string>& names,
                                 const std::vector<double>&      values,
                                 const std::vector<double>&      errors)
{
    _parameters = std::make_unique<ROOT::Minuit2::MnUserParameters>(values, errors);
    for (size_t i = 0; i < names.size(); i++) {
        _parameters->SetName(i, names[i]);
    }
}

std::vector<double> CharmFitterBase::ratios(void) const
{

    auto rsCounts{this->getRSBinContent()};
    auto wsCounts{this->getWSBinContent()};

    std::vector<double> ratio(_numBins);
    for (size_t i = 0; i < _numBins; ++i) {
        ratio[i] = wsCounts[i] / rsCounts[i];
    };

    return ratio;
}

std::vector<double> CharmFitterBase::errors(void) const
{
    auto rsCounts{this->getRSBinContent()};
    auto wsCounts{this->getWSBinContent()};
    auto ratio{this->ratios()};

    std::vector<double> error(_numBins);
    // Errors in a quotient add in quadrature; assuming our datasets obey Poisson statistics, we find that the
    // fractional error in our ratio is sqrt((N+M)/NM)
    for (size_t i = 0; i < _numBins; ++i) {
        error[i] = std::sqrt((wsCounts[i] + rsCounts[i]) / (wsCounts[i] * rsCounts[i])) * ratio[i];
    }

    return error;
}

TMatrixD CharmFitterBase::_covarianceVector2CorrelationMatrix(const std::vector<double>& covarianceVector,
                                                              const std::vector<double>& fitErrors)
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
    std::vector<double> errors = fitErrors;
    std::vector<size_t> paramsToRemove{};
    for (auto it = fixedParams.begin(); it != fixedParams.end(); ++it) {
        paramsToRemove.push_back(it->Number());
    }

    std::sort(paramsToRemove.rbegin(), paramsToRemove.rend());
    for (auto it = paramsToRemove.begin(); it != paramsToRemove.end(); ++it) {
        errors.erase(errors.begin() + *it);
    }

    // We need to divide each element in our vector of covariances with the standard deviation two parameters
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

FitResults_t CharmFitterBase::_performFit()
{
    // Check that we have no zero, inf or nan in our data
    for (const auto x : this->ratios()) {
        // Use TMath::Finite instead of std::is_finite() because for some reason including the ROOT headers makes
        // std::is_finite(-inf) sometimes return true...
        if (!TMath::Finite(x) || x == 0.0) {
            std::cerr << "Data must be finite: value " << x << " encountered." << std::endl;
            throw D2K3PiException();
        }
    }

    // Create a minimiser
    ROOT::Minuit2::MnMigrad migrad(*_fitFcn, *_parameters, 2U);

    // Minimuse chi squared as defined by our _fitFcn
    ROOT::Minuit2::FunctionMinimum min{migrad()};

    // Check that our solution is "valid"
    // I think this checks that the call limit wasn't reached and that the fit converged, though it's never possible
    // to be sure with Minuit2
    if (!min.IsValid()) {
        throw BadFitException(min);
    }

    // Fill in the fit parameters' values, errors, correlation matrix and fit statistic value.
    // The fit function will need to be set by the child class, as we don't know its form at this stage
    return FitResults_t{
        min.Fval(),
        min.UserParameters().Params(),
        min.UserParameters().Errors(),
        _covarianceVector2CorrelationMatrix(min.UserCovariance().Data(), min.UserParameters().Errors())};
}

std::vector<std::vector<double>> twoDParamScan(CharmFitter::CharmFitterBase&        Fitter,
                                               const std::function<double(double)>& efficiency,
                                               const std::string&                   param1,
                                               const std::string&                   param2,
                                               const std::vector<double>&           range1,
                                               const std::vector<double>&           range2)
{
    // Progress bar
    std::cout << "Performing " << param1 << "-" << param2 << " scan..." << std::endl;
    boost::progress_display progressBar(range1.size() * range2.size());

    // Fix the parameters we're interested in
    Fitter.fixParameter(param1);
    Fitter.fixParameter(param2);

    std::vector<std::vector<double>> result(range1.size(), std::vector<double>(range2.size(), 0.0));
    for (size_t i = 0; i < range1.size(); ++i) {
        Fitter.setParameter(param1, range1[i]);
        for (size_t j = 0; j < range2.size(); ++j) {
            Fitter.setParameter(param2, range2[j]);
            result[i][j] = Fitter.fit(efficiency).fitStatistic;

            ++progressBar;
        }
    }

    return result;
}

void savePlot(const CharmFitterBase&      Fitter,
              TF1&                        bestFitFunction,
              const std::string&          path,
              const std::string&          title,
              const util::LegendParams_t& legend)
{
    // Find bin errors (width/2)
    std::vector<double> binErrors = Fitter.getBinWidths();
    std::transform(binErrors.begin(), binErrors.end(), binErrors.begin(), [](const double x) { return x / 2; });

    size_t numBins = binErrors.size();

    // Create + save plot
    TGraphErrors graph(
        numBins, Fitter.getBinCentres().data(), Fitter.ratios().data(), binErrors.data(), Fitter.errors().data());

    graph.SetTitle(std::string{title + ";time/ns;DCS/CF ratio"}.c_str());
    util::saveObjectsToFile<TGraphErrors>(
        std::vector<TObject*>{&graph, &bestFitFunction}, {"AP", "CSAME"}, {"Data", "Best Fit"}, path, legend);
}

} // namespace CharmFitter
