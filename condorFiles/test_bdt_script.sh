#!/bin/bash
# Check that I can run my BDT script on Condor

set -e

LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH

source /software/mh19137/venvs/python3/bin/activate
which python
locate libc.so.6

echo "Running test BDT script"
python bdt_reweighting_test.py

