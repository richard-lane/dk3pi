sWeighting
==========

Our data is a combination of signal + background; we want to find a way to extract the signal.

sPlot method: input a signal + background model; fit to data; reweight data (sWeight)

 * Fit using discriminating variables {**y<sub>i</sub>**}
 * Project out some variable **x<sub>i</sub>** not in {**y<sub>i</sub>**}
 * Plot a histogram of all events in the dataset, each appropriately weighted
 * If the pdf for **x** is known then the sPlot can be checked to see if the fit has correctly determined it.


 Note: turns out it's against the fundamental idea of ROOT files to have my sWeighting code write some new branches in the existing file. What I've done instead is write a small ROOT file containing only the sWeighting branches, maybe these should be combined with hadd or something
