"""
Do a check AmpGen scan

"""
import numpy as np
import sys, os
import matplotlib.pyplot as plt
import common

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "build/charmFitter/python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "dk3pi-efficiency-model"))

import libcleoScan
from model.creation import utils

def main():
    # Read AmpGen
    rs_t = utils.read_ampgen("rs_Dbar02piKpipi_opt.root")[16]
    ws_t = utils.read_ampgen("ws_D02piKpipi_opt.root")[16]

    # Find out how many WS evts to keep
    x, y = libcleoScan.WORLD_AVERAGE_X, libcleoScan.WORLD_AVERAGE_Y
    r = 0.055
    re_z, im_z = 0.0, 0.0
    width = 2.5

    decay_params = [x, y, r, im_z, re_z, width]

    print(f"max RS:\t{np.max(rs_t)}\nmax WS:\t{np.max(ws_t)}")
    max_time = np.max(rs_t)

    num_ws = libcleoScan.numWS(len(rs_t), decay_params, max_time)
    ws_keep_fraction = num_ws / len(ws_t)
    ws_keep_mask = np.random.random(len(ws_t)) < ws_keep_fraction

    ws_t = ws_t[ws_keep_mask]
    print(f"Ideal num WS: {num_ws}\nActual: {len(ws_t)}")

    bins = np.logspace(np.log10(0.01), np.log10(max_time), 5)
    bins = np.insert(bins, 0, 0.0)
    print(bins)

    # Set up params
    n_re_z, n_im_z = 30, 30
    reZVals = np.linspace(-0.9, 0.9, n_re_z)
    imZVals = np.linspace(-0.9, 0.9, n_im_z)
    rsWeights    = np.ones(len(rs_t))
    wsWeights    = np.ones(len(ws_t))
    initialErrs  = [1, 1, 1, 1, 1, 1]

    # Create a fitter
    Fitter = libcleoScan.ConstrainedFitter(bins, decay_params, initialErrs)
    Fitter.addRSPoints(rs_t, np.ones_like(rs_t))
    Fitter.addWSPoints(ws_t, np.ones_like(ws_t))
    Fitter.fixParameters(("width", "z_im", "z_re"))

    print(Fitter.ratios())
    l = []
    for r in reZVals:
        for i in imZVals:
            Fitter.setParameter("z_re", r)
            Fitter.setParameter("z_im", i)
            l.append(Fitter.fit(lambda x: 1).fitStatistic)
            common.reset_params(Fitter, decay_params)
    l = np.array(l)

    # Subtract off minimum
    l = 2 * np.log(np.array(l))  # Should be -2?
    l -= np.min(l)
    l = np.sqrt(l)

    # Convert to 2d array
    l = l.reshape((n_re_z, n_im_z), order="C")

    # Plot something
    print(l)
    plt.contourf(reZVals, imZVals, l, [0, 1, 2, 3])
    plt.colorbar()

    path = "ampgen_scan.png"
    print(f"saving {path}")
    plt.savefig(path)


if __name__ == "__main__":
    main()

