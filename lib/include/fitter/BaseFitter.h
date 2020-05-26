#ifndef BASEFITTER_H
#define BASEFITTER_H

#include "FitterUtils.h"

#include "TGraphErrors.h"

/*
 * Base class that other fitters inherit from
 *
 * Holds the result of the fit as FitResults_t; also has attributes
 * for a plot and chi-squared/likelihood statistic.
 */
class BaseFitter
{
  public:
    /*
     * Parameters describing the fit
     */
    FitResults_t fitParams;

    /*
     * ROOT TGraph object, used for saving plots and stuff
     */
    std::unique_ptr<TGraphErrors> plot = nullptr;

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

#endif // BASEFITTER_H
