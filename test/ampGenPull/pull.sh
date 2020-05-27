#!/bin/bash

# Must have built target "ampgenpull" in the "build/test/ampGenPull" dir for this to work

set -e

DIR="$(dirname $(readlink -f $0))"

# Find how many of each event type to generate
# This is based on making some numbers up.
MEANNUMDCS=2313
MEANNUMCF=700000

# Generate ROOT files for D and Dbar events
# At the moment this generates equal numbers of both but this will need to change...
$AMPGENROOT/build/bin/Generator $DIR/../../AmpGenTools/options/Dbar02piKpipi.opt --nEvents $MEANNUMCF --EventType "Dbar0 K+ pi- pi- pi+" --Output dBar.root --GenerateTimeDependent
$AMPGENROOT/build/bin/Generator $DIR/../../AmpGenTools/options/generate_mixing.opt --nEvents $MEANNUMDCS --EventType "D K+ pi- pi- pi+" --Output d.root

# Read the ROOT files, calculate ratio of decay counts, perform a fit and save the fit parameters to a file
$DIR/../../build/test/ampGenPull/ampgenpull d.root dBar.root
