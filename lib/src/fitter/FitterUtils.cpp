#include <float.h>
#include <iostream>

#include "D2K3PiError.h"
#include "FitterUtils.h"

#include "TMath.h"

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