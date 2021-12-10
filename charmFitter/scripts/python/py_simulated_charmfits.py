"""
Do a check for charm fitter scan using the python binding

"""
# Set up include dir
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), "build/charmFitter/python"))

import cleoScan
import libcleoScan
import common


def _normalise(l: np.ndarray):
    """
    Arrays are passed by reference, effectively

    """
    # Subtract off minimum
    l[l != 0] -= np.min(l[l != 0])

    # Take sqrt to convert to sigmas
    return np.sqrt(l)


def main(phsp_bin):
    # What I expect Z should be, roughly, from the CLEO likelihood
    best_zs = {0: (-0.1, 0.7), 1: (-0.6, 0.6), 2: (-0.5, 0.2), 3: (0.25, 0.0)}

    # Define parameters used for simulating decays
    width = 2500.0
    decay_params = [0.0039183, 0.0065139, 0.055, best_zs[phsp_bin][1], best_zs[phsp_bin][0], width]

    max_time = 0.004
    n_rs = 10_000_000  # Number of RS evts to generate
    bins = [0, 0.00003, 0.00008, 0.00015, 0.00020, 0.00030, 0.00045, 0.00058, 0.00060, 0.00085, 0.0012, 0.004]

    # Simulate decays
    seed = np.random.randint(0, np.iinfo(np.uint32).max)
    rs_t, ws_t = libcleoScan.simulate(n_rs, decay_params, max_time, seed)

    # Set up params
    n_re_z, n_im_z = 50, 50
    reZVals = np.linspace(-0.99, 0.99, n_re_z)
    imZVals = np.linspace(-0.99, 0.99, n_im_z)
    rsWeights = np.ones(len(rs_t))
    wsWeights    = np.ones(len(ws_t))
    initialErrs  = [1, 1, 1, 1, 1, 1]

    # Find chisqs
    charm_likelihood = _normalise(cleoScan.charmScan(reZVals,
                                          imZVals,
                                          rs_t,
                                          rsWeights,
                                          ws_t,
                                          wsWeights,
                                          bins,
                                          decay_params,
                                          initialErrs))

    cleo_likelihood = cleoScan.cleoLikelihoods(reZVals, imZVals, decay_params, phsp_bin)

    # Find the min non-NaN value, replace NaNs with it
    cleo_max = np.min(cleo_likelihood[~np.isnan(cleo_likelihood)])  # Min instead of max because we're dealing with -2LL
    cleo_likelihood[np.isnan(cleo_likelihood)] = cleo_max
    cleo_likelihood = _normalise(-2.0 * cleo_likelihood)

    combined_likelihood = _normalise(cleoScan.cleoScan(reZVals,
                                            imZVals,
                                            rs_t,
                                            rsWeights,
                                            ws_t,
                                            wsWeights,
                                            bins,
                                            decay_params,
                                            initialErrs,
                                            phsp_bin))

    # Plot something
    titles = "Charm Mixing", "CLEO", "Both"
    fig, ax = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
    for a, l, t in zip(ax, (charm_likelihood, cleo_likelihood, combined_likelihood), titles):
        contours = a.contourf(reZVals, imZVals, l, list(range(4)), cmap="turbo")
        a.add_patch(Circle((0, 0), 1.0, color="k", linestyle="--", fill=False))
        a.set_title(t)
        a.set_xlabel(r"Re(Z)")
        a.set_ylabel(r"Im(Z)")
        a.set_xlim(-1.0, 1.0)
        a.set_ylim(-1.0, 1.0)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(contours, cax=cbar_ax)
    cbar.ax.set_ylabel(r"$\sigma$", rotation=90)

    fig.suptitle("$Phsp bin" + str(phsp_bin) + " : Z_\Omega^{K3\pi} Constraint$")

    # Expected relationship between components of Z
    decay_params = [0.0039183, 0.0065139, 0.055, 0.5, 0.0, width]
    gradient = -decay_params[1] / decay_params[0]
    intercept = decay_params[3] - gradient * decay_params[4]
    expected_im_z = [intercept + gradient*x for x in reZVals]
    for a in (ax[0], ax[2]):
        a.plot(reZVals, expected_im_z, "k--")

    plt.savefig(f"py_scans_bin{phsp_bin}.png")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        phsp_bin = 0
    else:
        phsp_bin = int(sys.argv[1])

    main(phsp_bin)

