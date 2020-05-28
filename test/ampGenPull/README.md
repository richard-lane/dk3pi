Pull study using decays generated with AmpGen
This might take ages to run, depending on the number of events used

Python script creates events with AmpGen; C++ file performs the time-binning, takes ratios, performs fit + saves data to a text file.

Must build target ampgenpull in build/test/ampGenPull/ before running
the pull study python script.
