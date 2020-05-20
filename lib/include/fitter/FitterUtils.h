/*
 * Utility structs and stuff for fitting
 */
#ifndef FITTERUTILS_H
#define FITTERUTILS_H

#include <vector>

#include "TMatrixD.h"

/*
 * The data needed to perform the fit
 * Bin centres, bin widths, our data and its error
 *
 * Parameters are all initialised to zero by default and should be set by the user.
 */
typedef struct FitData {

    /*
     * Default constructor sets all attributes to zero
     */
    FitData(void);

    /*
     * Check that all parameters are vectors of the same length, our bins don't overlap and data is all finite.
     * Note that this will not accept a dataset which contains zeros.
     * Sets numPoints
     */
    FitData(const std::vector<double>& myBinCentres,
            const std::vector<double>& myBinWidths,
            const std::vector<double>& myData,
            const std::vector<double>& myErrors);

    std::vector<double> binCentres{0.0};

    /*
     * This is not the total width of the bins; it is the error in the bin, or half the width
     */
    std::vector<double> binErrors{0.0};

    std::vector<double> data{0.0};
    std::vector<double> errors{0.0};

    /*
     * The number of data points that our dataset represents
     * Set by FitData constructor.
     */
    size_t numPoints{0};
} FitData_t;

/*
 * The parameters describing how a function fits data
 */
typedef struct PolynomialFitResults {

    /*
     * Vector of the fit parameters
     */
    std::vector<double> fitParams{};

    /*
     * Vector of the errors in our fit parameters
     */
    std::vector<double> fitParamErrors{};

    /*
     * Pointer to correlation matrix of fit parameters
     */
    std::unique_ptr<TMatrixD> correlationMatrix{};
} FitResults_t;

#endif // FITTERUTILS_H
