AmpGenTools
===========

This dir is for extra tools + stuff pertaining to AmpGen, such as options files (*.opt)
It also contains the D-> K3pi CF and DCS amplitude model file things; these are needed to perform the phase space binning (lib/k3pi_binning.h).

How to use options files to generate D->K3pi
--------------------------------------------

(We're looking for D -> K+ pi- pi- pi+ decays)
We expect the D0 amplitude to have contributions from both the DCS direct path and the mixing then CF decay path, so this amplitude depends
quite strongly on mixing.
We expect the Dbar0 amplitude to be completely dominated by the CF direct path, so this one doesn't depend on mixing.

The AmpGen .opt files generate_mixing.opt and Dbar02piKpipi.opt model these two situations respectively.

To generate D0 events (with mixing), use:
> $AMPGENROOT/build/bin/Generator AmpGenTools/options/generate_mixing.opt --nEvents 1000 --EventType "D0 K+ pi- pi- pi+"
> --Output d.root

To generate D events (no mixing), use:
> $AMPGENROOT/build/bin/Generator AmpGenTools/options/Dbar02piKpipi.opt --nEvents 1000 --EventType "Dbar0 K+ pi- pi- pi+"
> --Output dBar.root --GenerateTimeDependent

Obviously change the number of events and filename to whatever you want

