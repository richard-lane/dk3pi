#include <float.h>
#include <iostream>

#include "D2K3PiError.h"
#include "FitterUtils.h"

#include "TMath.h"

FitData::FitData()
{
    ;
}

FitData::FitData(const std::vector<double>& myBinLimits,
                 const std::vector<double>& myData,
                 const std::vector<double>& myErrors)
{

    // Check that our bin limitsare sorted
    if (!std::is_sorted(myBinLimits.begin(), myBinLimits.end())) {
        std::cerr << "Bins should be sorted" << std::endl;
        throw D2K3PiException();
    }

    // Check that all of our vectors are the same length
    numBins           = myBinLimits.size() - 1;
    size_t dataLength = myData.size();
    size_t errorsSize = myErrors.size();
    if (numBins != errorsSize || dataLength != errorsSize) {
        std::cerr << "FitData parameters must all be the same length." << std::endl;
        throw D2K3PiException();
    }

    // Check that the data is all finite because -infs and 0s keep sneaking in and its annoying
    for (auto it = myData.begin(); it != myData.end(); ++it) {
        // Use TMath::Finite instead of std::is_finite() because for some reason including the ROOT headers makes
        // std::is_finite(-inf) sometimes return true...
        if (!TMath::Finite(*it) || *it == 0.0) {
            std::cerr << "Data must be finite: value " << *it << " encountered." << std::endl;
            throw D2K3PiException();
        }
    }

    // Now that all the checks have passed, set the class attributes
    binLimits = myBinLimits;
    data      = myData;
    errors    = myErrors;

    binCentres = std::vector<double>(numBins);
    for (size_t i = 0; i < numBins; ++i) {
        binCentres[i] = (binLimits[i] + binLimits[i + 1]) / 2;
    }
}

MyBaseFcn::MyBaseFcn(const std::vector<double>& data,
                     const std::vector<double>& times,
                     const std::vector<double>& errors,
                     const IntegralOptions_t&   integralOptions)
    : theMeasurements(data), thePositions(times), theMVariances(errors), theErrorDef(1.),
      _integralOptions(integralOptions)
{
    ;
}

MyBaseFcn::~MyBaseFcn()
{
    ;
}

double MyBaseFcn::Up() const
{
    return theErrorDef;
}

std::vector<double> MyBaseFcn::measurements() const
{
    return theMeasurements;
}

std::vector<double> MyBaseFcn::positions() const
{
    return thePositions;
}

std::vector<double> MyBaseFcn::variances() const
{
    return theMVariances;
}

void MyBaseFcn::setErrorDef(double def)
{
    theErrorDef = def;
}
