"""
Do a check AmpGen scan

"""
import numpy as np
import sys, os
import matplotlib.pyplot as plt
import common
import pickle

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), "build/charmFitter/python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))), "dk3pi-efficiency-model"))

import libcleoScan
from model.creation import utils
from model.util import definitions


def main():
    # Read AmpGen
    rs_t = utils.read_ampgen("rs_Dbar02piKpipi_opt.root")[16]
    ws_t = utils.read_ampgen("ws_D02piKpipi_opt.root")[16]

    # Read efficiency weights
    with open("rs_weights.pickle", "rb") as f:
        rs_weights = pickle.load(f)
    with open("ws_weights.pickle", "rb") as f:
        ws_weights = pickle.load(f)

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
    ws_weights = ws_weights[ws_keep_mask]
    print(f"Ideal num WS: {num_ws}\nActual: {len(ws_t)}")

    ws_t = np.delete(ws_t, np.where(ws_weights == 0.0))
    ws_weights = ws_weights[ws_weights != 0.0]

    rs_t = rs_t[rs_t > definitions.MIN_TIME_PS]

    print(f"After deleting: {len(ws_t)}")

    bins = np.logspace(np.log10(definitions.MIN_TIME_PS), np.log10(max_time), 5)
    print(bins)

    # Set up params
    n_re_z, n_im_z = 30, 30
    reZVals = np.linspace(-0.9, 0.9, n_re_z)
    imZVals = np.linspace(-0.9, 0.9, n_im_z)
    initialErrs  = [1, 1, 1, 1, 1, 1]

    # Create a fitter
    UnweightedFitter = libcleoScan.ConstrainedFitter(bins, decay_params, initialErrs)
    UnweightedFitter.addRSPoints(rs_t, np.ones_like(rs_t))
    UnweightedFitter.addWSPoints(ws_t, np.ones_like(ws_t))
    UnweightedFitter.fixParameters(("width", "z_im", "z_re"))

    WeightedFitter = libcleoScan.ConstrainedFitter(bins, decay_params, initialErrs)
    WeightedFitter.addRSPoints(rs_t, rs_weights)
    WeightedFitter.addWSPoints(ws_t, ws_weights)
    WeightedFitter.fixParameters(("width", "z_im", "z_re"))

    print(f"{UnweightedFitter.ratios()=}")
    print(f"{WeightedFitter.ratios()=}")

    weighted_likelihood = []
    unweighted_likelihood = []
    for r in reZVals:
        for i in imZVals:
            UnweightedFitter.setParameter("z_re", r)
            UnweightedFitter.setParameter("z_im", i)
            unweighted_likelihood.append(UnweightedFitter.fit(lambda x: 1).fitStatistic)

            WeightedFitter.setParameter("z_re", r)
            WeightedFitter.setParameter("z_im", i)
            weighted_likelihood.append(WeightedFitter.fit(lambda x: 1).fitStatistic)

            common.reset_params(UnweightedFitter, decay_params)
            common.reset_params(WeightedFitter, decay_params)

    weighted_likelihood = np.array(weighted_likelihood)
    unweighted_likelihood = np.array(unweighted_likelihood)

    # Subtract off minimum
    weighted_likelihood = 2 * np.log(np.array(weighted_likelihood))  # Should be -2?
    weighted_likelihood -= np.min(weighted_likelihood)

    unweighted_likelihood = 2 * np.log(np.array(unweighted_likelihood))  # Should be -2?
    unweighted_likelihood -= np.min(unweighted_likelihood)

    # Convert to 2d array
    weighted_likelihood = weighted_likelihood.reshape((n_re_z, n_im_z), order="C")
    unweighted_likelihood = unweighted_likelihood.reshape((n_re_z, n_im_z), order="C")

    # Plot something
    print(unweighted_likelihood)
    print(weighted_likelihood)

    fig, ax = plt.subplots(1, 2)
    ax[0].contourf(reZVals, imZVals, unweighted_likelihood, [0, 1, 2, 3])
    ax[1].contourf(reZVals, imZVals, weighted_likelihood, [0, 1, 2, 3])

    ax[0].set_title("Unweighted")
    ax[1].set_title("Weighted")

    path = "ampgen_scan.png"
    print(f"saving {path}")
    plt.savefig(path)


if __name__ == "__main__":
    main()

