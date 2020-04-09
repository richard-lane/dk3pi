/*
 * Base class that other fitters inherit from
 *
 * Holds the result of the fit as FitResults_t; also has attributes
 * for a plot and chi-squared/likelihood statistic.
 */

#include <FitterUtils.h>

#include "TGraphErrors.h"

class BaseFitter
{
  public:
    /*
     * Parameters describing the fit
     */
    FitResults_t fitParams;

    /*
     * ROOT TGraph object used for holding input data.
     */
    std::unique_ptr<TGraphErrors> plot = nullptr;

    /*
     * The test statistic that is optimised when running our test, e.g. chi squared or a likelihood
     */
    std::unique_ptr<double> statistic = nullptr;

  protected:
    /*
     * The data to be fit.
     */
    BaseFitter(const FitData_t& fitData);

    /*
     * The data used to make the fit
     */
    FitData_t _fitData;
};