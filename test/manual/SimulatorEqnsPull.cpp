/*
#include "MinuitFcns.h"
#include "../pull_study/PullStudyHelpers.h"
#include "physics.h"
#include "util.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"

#include <boost/progress.hpp>


 * Fitter class for fitting to (a + bt + ct^2)e^-(width * t)
 *
 * Performs an unbinned likelihood fit using Minuit2 API
/
class TestSimulatorFitter : public BasePolynomialFcn
{
  public:
    TestSimulatorFitter(const std::vector<double> &decayTimes, const double width, const double maxTime)
        : BasePolynomialFcn(decayTimes, std::vector<double>{}, std::vector<double>{}), _width(width), _maxTime(maxTime)
    {
        ;
    }

    ~TestSimulatorFitter() { ; }

    virtual double operator()(const std::vector<double> &parameters) const
    {

        double integral      = Phys::analyticalDcsIntegral(0, _maxTime, parameters, _width);
        double logLikelihood = 0.0;
        for (size_t i = 0; i < theMeasurements.size(); ++i) {
            double prob = (parameters[0] + parameters[1] * theMeasurements[i] +
                           parameters[2] * theMeasurements[i] * theMeasurements[i]) *
                          std::exp(-_width * theMeasurements[i]) / integral;
            logLikelihood -= 2 * std::log(prob);
        }
        return logLikelihood;
    }

  private:
    double _width{0.0};
    double _maxTime{0.0};
};

*
 * Overlay a plot of a fit onto data
 *
void plotFit(ROOT::Minuit2::FunctionMinimum min,
             SimulatedDecays                MyDecays,
             const double                   maxTime,
             const double                   width,
             const size_t                   numEvents)
{
    double integral = Phys::analyticalDcsIntegral(0, maxTime, min.UserParameters().Params(), width);
    auto prob = [&](double t) { return Phys::wrongSignDecayRate(t, min.UserParameters().Params(), width) / integral; };

    // Bins
    size_t              numBins  = 50;
    double              binWidth = maxTime / numBins;
    std::vector<double> binCentres(numBins);
    std::vector<double> binLimits(numBins);
    for (size_t i = 0; i <= numBins; ++i) {
        binLimits[i] = i * maxTime / numBins;
    }

    TH1D *hist = new TH1D("fit", "fit", numBins, binLimits.data());
    hist->FillN(MyDecays.WSDecayTimes.size(), MyDecays.WSDecayTimes.data(), nullptr);

    // Model values in bin centres
    std::vector<double> model(numBins);
    std::vector<double> error(numBins);
    for (size_t i = 0; i < numBins; ++i) {
        binCentres[i] = (binLimits[i] + binLimits[i + 1]) / 2;
        model[i]      = numEvents * prob(binCentres[i]) * binWidth;
        error[i]      = std::sqrt(hist->GetBinContent(i + 1));
    }

    // Create graph
    TGraphErrors *graph = new TGraphErrors(numBins, binCentres.data(), model.data(), 0, error.data());

    // Plot
    TCanvas *c = new TCanvas();
    c->cd();
    hist->Draw();
    graph->Draw("CSAME");
    c->SaveAs("fit.pdf");

    delete c;
    delete hist;
    delete graph;
}

void scan(const double               mean,
          const double               std,
          const std::vector<double> &expectedParams,
          const std::string &        fix,
          const double               width,
          const double               maxTime,
          const SimulatedDecays &    MyDecays)
{
    // Parameter scan
    size_t              n = 100;
    std::vector<double> aVals(n);
    std::vector<double> likelihoodVals(n);
    double              min = mean - 2 * std;
    double              max = mean + 2 * std;

    for (size_t i = 0; i < n; ++i) {
        // Scan +-2 sigma
        aVals[i] = min + ((double)i / (double)n) * (max - min);

        ROOT::Minuit2::MnUserParameters scanPar;
        scanPar.Add("a", expectedParams[0], 0.1);
        scanPar.Add("b", expectedParams[1], 0.1);
        scanPar.Add("c", expectedParams[2], 0.1);

        scanPar.SetValue(fix, aVals[i]);

        TestSimulatorFitter     scanFitter(MyDecays.WSDecayTimes, width, maxTime);
        ROOT::Minuit2::MnMigrad scanMigrad(scanFitter, scanPar);
        scanMigrad.Fix(fix.c_str());
        ROOT::Minuit2::FunctionMinimum scanMin = scanMigrad();
        likelihoodVals[i]                      = scanMin.Fval();
        std::cout << scanMin;
    }
    TGraph *Graph = new TGraph(n, aVals.data(), likelihoodVals.data());
    util::saveObjectToFile(Graph, (fix + "scan.pdf").c_str(), "AP");
    delete Graph;
}

 * Generate DCS and CF events, fit to them using (a + bt + ct^2)e^(-width*t)
void dcsTimesPull()
{
    DecayParams_t MyParams{
        .x     = 0.0037,
        .y     = 0.0066,
        .r     = 0.055,
        .z_im  = -0.2956,
        .z_re  = 0.7609,
        .width = 2439.0,
    };
    double maxTime = 5 / MyParams.width;

    // Create our Decay simulator object and generate CF and DCS times using accept-reject
    size_t numDecays      = 1e5;
    size_t numExperiments = 100;

    // Initialise expected parameters
    std::vector<double>             expectedParams = util::expectedParams(MyParams);
    ROOT::Minuit2::MnUserParameters dcsPar;
    dcsPar.Add("a", expectedParams[0], 0.1);
    dcsPar.Add("b", expectedParams[1], 0.1);
    dcsPar.Add("c", expectedParams[2], 0.1);

    std::vector<double> aPull(numExperiments, -1);
    std::vector<double> bPull(numExperiments, -1);
    std::vector<double> cPull(numExperiments, -1);

    boost::progress_display showProgress(numExperiments);

    for (size_t i = 0; i < numExperiments; ++i) {
        SimulatedDecays MyDecays = SimulatedDecays(maxTime, MyParams, 0);
        MyDecays.findDcsDecayTimes(numDecays);

        // Create fitter
        TestSimulatorFitter     dcsFitter(MyDecays.WSDecayTimes, MyParams.width, maxTime);
        ROOT::Minuit2::MnMigrad dcsMigrad(dcsFitter, dcsPar);
        dcsMigrad.Fix("a");

        // Perform fit
        ROOT::Minuit2::FunctionMinimum dcsMin = dcsMigrad();
        if (!dcsMin.IsValid()) {
            std::cout << dcsMin;
        }

        // Pull
        aPull[i] = (dcsMin.UserParameters().Params()[0] - expectedParams[0]) / dcsMin.UserParameters().Errors()[0];
        bPull[i] = (dcsMin.UserParameters().Params()[1] - expectedParams[1]) / dcsMin.UserParameters().Errors()[1];
        cPull[i] = (dcsMin.UserParameters().Params()[2] - expectedParams[2]) / dcsMin.UserParameters().Errors()[2];

        ++showProgress;

        // Plot
        // plotFit(dcsMin, MyDecays, maxTime, MyParams.width, numDecays);

        // scan(dcsMin.UserParameters().Params()[0],
        //      dcsMin.UserParameters().Errors()[0],
        //      expectedParams,
        //      "a",
        //      MyParams.width,
        //      maxTime,
        //      MyDecays);

        // scan(dcsMin.UserParameters().Params()[1],
        //      dcsMin.UserParameters().Errors()[1],
        //      expectedParams,
        //      "b",
        //      MyParams.width,
        //      maxTime,
        //      MyDecays);

        // scan(dcsMin.UserParameters().Params()[2],
        //      dcsMin.UserParameters().Errors()[2],
        //      expectedParams,
        //      "c",
        //      MyParams.width,
        //      maxTime,
        //      MyDecays);
    }

    PullStudyHelpers::plot_parameter_distribution("a", aPull, numExperiments);
    PullStudyHelpers::plot_parameter_distribution("b", bPull, numExperiments);
    PullStudyHelpers::plot_parameter_distribution("c", cPull, numExperiments);

    std::cout << PullStudyHelpers::meanAndStdDev(aPull).first << "+-" << PullStudyHelpers::meanAndStdDev(aPull).second
              << std::endl;
    std::cout << PullStudyHelpers::meanAndStdDev(bPull).first << "+-" << PullStudyHelpers::meanAndStdDev(bPull).second
              << std::endl;
    std::cout << PullStudyHelpers::meanAndStdDev(cPull).first << "+-" << PullStudyHelpers::meanAndStdDev(cPull).second
              << std::endl;
}
*/

int main()
{
    //dcsTimesPull();
    return 0;
}