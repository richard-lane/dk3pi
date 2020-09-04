Efficiency
====

Particle detectors aren't perfect- each event has a certain probability
of being detected which is in general less than 1.
For this analysis (D -> K pi pi pi) we have 5-dimensional final state
phase space, so our efficiency will be a function over a 5d space.

We can estimate our efficiency from Monte Carlo data- it (presumably) reproduces
the detector quite well. However, it won't be possible to generate enough data
to find the efficiency directly from Monte Carlo- if we have 100 bins along each of
our 5 axes, that will be 10^10 bins which is far too many to find an efficiency value
in each one from data.

Instead, we have to do something clever...

 - [Chow Liu](./chowLiu/)
