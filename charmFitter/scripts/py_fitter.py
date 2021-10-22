"""
Do a check for CLEO scan using the python binding

"""
# Set up include dir
import numpy as np
import sys, os
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "dk3pi-efficiency-model"))

import cleoScan
from model.creation import utils

def main():
    # Read AmpGen
    rs_t = utils.read_ampgen("rs_Dbar02piKpipi_opt.root")[16]
    ws_t = utils.read_ampgen("ws_D02piKpipi_opt.root")[16]
    ws_t = ws_t[0:200000]


    bins = np.linspace(0, 5, 25)
    bins = np.append(bins, 6.0)

    # kw = {"bins": bins, "histtype": "step"}
    # plt.hist(rs_t, **kw)
    # plt.hist(ws_t, **kw)
    # plt.yscale("log")
    # plt.savefig("tmp.png")

    # Set up params
    reZVals = np.linspace(-0.9, 0.9, 30)
    imZVals = np.linspace(-0.9, 0.9, 30)
    rsWeights    = np.ones(len(rs_t))
    wsWeights    = np.ones(len(ws_t))
    initialVals  = [1, 0, 0, 0, 0, 0.41]
    initialErrs  = [1, 1, 1, 1, 1, 1]
    binNumber    = 0

    # Call fcn
    l = cleoScan.cleoScan(reZVals,
                          imZVals,
                          rs_t,
                          rsWeights,
                          ws_t,
                          wsWeights,
                          bins,
                          initialVals,
                          initialErrs,
                          binNumber)

    # Convert to 2d array
    l = l.reshape((len(reZVals), len(imZVals)), order="C")
    print(l)

    # Plot something
    plt.imshow(l)
    plt.savefig("tmp.png")


if __name__ == "__main__":
    main()

