"""
Do a fit to the AmpGen modeel with charm mixing turned on (generate_mixing.opt)

"""
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import common
import pickle

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), "build/charmFitter/python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))), "dk3pi-efficiency-model"))

import libcleoScan
from model.creation import utils
from model.util import definitions
from model.util import phsp_binning


def main():
    # Read AmpGen
    rs = utils.read_ampgen("rs_Dbar02piKpipi_opt.root")
    ws = utils.read_ampgen("mixed_generate_mixing_opt.root")

    rs_k, rs_pi1, rs_pi2, rs_pi3, rs_t = rs[0:4], rs[4:8], rs[8:12], rs[12:16], rs[16]
    ws_k, ws_pi1, ws_pi2, ws_pi3, ws_t = ws[0:4], ws[4:8], ws[8:12], ws[12:16], ws[16]

    # Find bins
    rs_bins = phsp_binning.bin_numbers(rs_k, rs_pi1, rs_pi2, rs_pi3, +1)
    ws_bins = phsp_binning.bin_numbers(ws_k, ws_pi1, ws_pi2, ws_pi3, +1)

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
    ws_bins = ws_bins[ws_keep_mask]
    print(f"Ideal num WS: {num_ws}\nActual: {len(ws_t)}")

    # Throw away low decay times...
    rs_bins = rs_bins[rs_t > definitions.MIN_TIME_PS]
    ws_bins = ws_bins[ws_t > definitions.MIN_TIME_PS]
    rs_t = rs_t[rs_t > definitions.MIN_TIME_PS]
    ws_t = ws_t[ws_t > definitions.MIN_TIME_PS]
    print(f"After veto: RS: {len(rs_t)}\nWS: {len(ws_t)}")
    bins = np.array([definitions.MIN_TIME_PS, 2.0, 3.0, 5.8])
    path = "mixed_ampgen_scan_veto.png"

    # Don't throw away low decay times
    # bins = np.array([0.0, definitions.MIN_TIME_PS, 2.0, 3.0, 5.8])
    # path = "mixed_ampgen_scan.png"

    # Set up params
    n_re_z, n_im_z = 30, 30
    reZVals = np.linspace(-0.9, 0.9, n_re_z)
    imZVals = np.linspace(-0.9, 0.9, n_im_z)
    initialErrs  = [1, 1, 1, 1, 1, 1]

    # Create a fitter
    Fitter = libcleoScan.ConstrainedFitter(bins, decay_params, initialErrs)
    Fitter.addRSPoints(rs_t[rs_bins == 1], np.ones_like(rs_t[rs_bins == 1]))
    Fitter.addWSPoints(ws_t[ws_bins == 1], np.ones_like(ws_t[ws_bins == 1]))
    Fitter.fixParameters(("width", "z_im", "z_re"))

    print(f"{Fitter.ratios()=}")

    likelihood = []
    for r in reZVals:
        for i in imZVals:
            Fitter.setParameter("z_re", r)
            Fitter.setParameter("z_im", i)
            likelihood.append(Fitter.fit(lambda x: 1).fitStatistic)

            common.reset_params(Fitter, decay_params)

    likelihood = np.array(likelihood)

    # Subtract off minimum
    likelihood = 2 * np.log(np.array(likelihood))  # Should be -2?
    likelihood -= np.min(likelihood)

    # Convert to 2d array
    likelihood = likelihood.reshape((n_re_z, n_im_z), order="C")
    print(np.min(likelihood))

    # Plot something
    fig, ax = plt.subplots()
    # ax.contourf(reZVals, imZVals, likelihood, [0, 1, 2, 3])
    img = ax.imshow(likelihood, vmin=0.0, vmax=3.0)
    cax = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.05)
    fig.colorbar(img, cax=cax, orientation="vertical")

    ax.set_yticks([])
    ax.set_xticks([])

    ax.set_title(f"AmpGen Mixed")

    print(f"saving {path}")
    plt.savefig(path)


if __name__ == "__main__":
    main()

