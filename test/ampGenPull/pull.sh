#!/bin/bash

set -e

DIR="$(dirname $(readlink -f $0))"

# Generate ROOT files for D and Dbar events
# At the moment this generates equal numbers of both but this will need to change...
$AMPGENROOT/build/bin/Generator $DIR/../../AmpGenTools/options/Dbar02piKpipi.opt --nEvents 1000 --EventType "Dbar0 K+ pi- pi- pi+" --Output dBar.root --GenerateTimeDependent
$AMPGENROOT/build/bin/Generator AmpGenTools/options/generate_mixing.opt --nEvents 1000 --EventType "D0 K+ pi- pi- pi+" --Output d.root

# Read the ROOT files, calculate ratio of decay counts, perform a fit and save the fit parameters to a file

