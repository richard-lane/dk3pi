#!/bin/bash
# This should probably be a python script

# Must have built target "ampgenpull" in the "build/test/ampGenPull" dir for this to work

set -e

DIR="$(dirname $(readlink -f $0))"

# This is based on making some numbers up
# Event numbers to use get written to a file
export MEANNUMDCS=231
export MEANNUMCF=70000

for i in {1..2}
do
  # Find how many of each event type to generate
  python $DIR/find_event_nums.py
  NUMDCS=$(sed '1q;d' event_nums.txt)
  NUMCF=$(sed '2q;d' event_nums.txt)

  # Generate ROOT files for D and Dbar events
  # At the moment this generates equal numbers of both but this will need to change...
  $AMPGENROOT/build/bin/Generator $DIR/../../AmpGenTools/options/Dbar02piKpipi.opt --nEvents $NUMCF --EventType "Dbar0 K+ pi- pi- pi+" --Output dBar.root --GenerateTimeDependent
  $AMPGENROOT/build/bin/Generator $DIR/../../AmpGenTools/options/generate_mixing.opt --nEvents $NUMDCS --EventType "D K+ pi- pi- pi+" --Output d.root

  # Read the ROOT files, calculate ratio of decay counts, perform a fit and save the fit parameters to a file
  $DIR/../../build/test/ampGenPull/ampgenpull d.root dBar.root
done

# Find the mean + std dev of pulls
python $DIR/find_pull.py
