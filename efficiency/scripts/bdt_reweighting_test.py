"""
BDT reweighting in pure python

"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# This is really horrible but i can't think of a better way of doing it
# Ideally i'd like to set a global python include path via the CMake build system...
sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
import reweighting
import reweight_utils


def main():
    # Find phsp points for prompt + SL datasets

    # Find weights for prompt + SL datasets

    # Split data into training + test data

    # Train the BDT on training data

    # Reweight the test prompt data

    # Compare the reweighted test prompt + test SL data by plotting some histograms
    pass


if __name__ == "__main__":
    main()
