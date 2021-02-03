# Multidimensional Goodness of Fit
In 1d, there are many methods for comparing distributions in 1d (\chi^2, KS test, ...)
Following on from [this](https://indico.cern.ch/event/975074/) meeting, a method for finding the goodness of fit of two multidimensional distributions is discussed here.

# Key Ideas
 - There are lots of way s of comparing distributions in 1d
 - There are comparatively fewer ways of comparing higher-dimensional distributions
 - ML tools give us ways of reducing a high-dimensional problem to a 1d one
 - Once we have reduced the dimensionality of our problem, we can run a test on the 1d distribution returned by the ML algorithm.
 - This also allows us to calculate a p-value for this 1d test.

 # Concrete example
  - Consider the phase space distribution for K->3\pi decays (5 dimensional)
  - We can simulate 5d distributions using an amplitude model ("model") or by using the full LHCb Monte Carlo ("MC")
  - How well do these distributions match in 5d?
  - First train a binary classifier to distinguish between model and MC events
  - Then use this classifier on model and MC holdout test sets to find the probability of categorising each as MC
  - Plot these probabilities on a histogram:
  - Find the \chi^2 score between these histograms and its corresponding p value
