"""
Do a check for charm fitter scan using the python binding

"""
# Set up include dir
import numpy as np
import sys, os
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), "build/charmFitter/python"))

import cleoScan
import libcleoScan
import common


def main():
    # Simulate decays
    # Define our parameters
    width = 2500.0
    decay_params = [0.0039183, 0.0065139, 0.055, 0.5, 0.0, width]

    max_time = 0.004
    n_rs = 1_000_000  # Number of RS evts to generate each time
    bins = [0, 0.00003, 0.00008, 0.00015, 0.00020, 0.00030, 0.00045, 0.00058, 0.00060, 0.00085, 0.0012, 0.004]

    # Simulate decays
    seed = np.random.randint(0, np.iinfo(np.uint32).max)
    rs_t, ws_t = libcleoScan.simulate(n_rs, decay_params, max_time, seed)

    # Set up params
    n_re_z, n_im_z = 25, 25
    reZVals = np.linspace(-1.0, 1.0, n_re_z)
    imZVals = np.linspace(-1.0, 1.0, n_im_z)
    rsWeights    = np.ones(len(rs_t))
    wsWeights    = np.ones(len(ws_t))
    initialErrs  = [1, 1, 1, 1, 1, 1]

    # Find chisqs
    l = cleoScan.cleoScan(reZVals,
             imZVals,
             rs_t,
             rsWeights,
             ws_t,
             wsWeights,
             bins,
             decay_params,
             initialErrs,
             0)

    # Subtract off minimum
    l[l != 0] -= np.min(l[l != 0])

    # Take sqrt to convert to sigmas
    l = np.sqrt(l)

    # Plot something
    plt.contourf(reZVals, imZVals, l, 10)
    plt.colorbar()
    plt.savefig("py_scan.png")


if __name__ == "__main__":
    main()

