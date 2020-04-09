#ifndef ROOTFITTER_H
#define ROOTFITTER_H

#include "BaseFitter.h"

/*
 * Fit to a polynomial (a + bt + ct^2) using ROOT's builtin TGraph Fit() method.
 */
class RootFitter : public BaseFitter
{
  public:
    /*
     * Calls parent constructor
     */
    RootFitter(const FitData_t& fitData);

    /*
     * Fit our data to the equation we expect to see, using ROOT's builtin TGraph fitter.
     * At the moment this is just a second order polynomial.
     *
     * minTime and maxTime define the range over which the function is defined.
     *
     * Populates fitParams and allocated memory to plot
     */
    void fit(const double minTime, const double maxTime, const std::string& options = "");

    /*
     * Save a plot of our data and best fit curve to file
     */
    void saveFitPlot(const std::string& plotTitle, const std::string& path);
};

#endif // ROOTFITTER_H
