#!/bin/bash
# Check that I can run my BDT script on Condor

echo "Running test BDT script"
which python
python --version
source /storage/mh19137/python_venv/bin/activate

which python
python --version
python bdt_reweighting_test.py

