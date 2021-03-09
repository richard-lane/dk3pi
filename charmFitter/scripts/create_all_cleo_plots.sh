#!/bin/bash

set -e

# Build the executable before executing

for BIN in {0..3}
do
  # Create plots
  ./charmFitter/scripts/cleo_combination.exe $BIN

  # Stitch together
  convert bin${BIN}_cleo.png bin${BIN}_mixing.png bin${BIN}_combined.png +append bin${BIN}_all.png

done
