Massively overcomplicated fitter

Uses Minuit2 APIs to perfom fits

Class structure:
 * BaseFitter
 * MinuitFitterBase
 * MinuitScannerBase

Are base clases that inherit from each other (top -> bottom)
They're a base class, base class for minuit fitting, base class for performing scans with Minuit respectively.

Then the actual user-facing fitters are:
 * RootFitter
 * MinuitPolynomialFitter
 * PhysicalFitter

For fit using ROOT, fit to (a + bt + ct2) using Minuit and fit to the eqn containing rD, x, y, Re(Z), Im(Z) etc. respectively.
They inherit from basefitter, minuitscannerbase and minuitscannerbase respectively.
