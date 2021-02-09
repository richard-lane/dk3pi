"""
Plot flat, RS and WS efficiency projections

"""
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
import plotting
import script_util
import reweight_utils


def _density_scale_factor(bins, data):
    """
    Find the scale factor used when plotting a plt hist with the density option

    """
    orig, _, _ = np.histogram(data, bins=bins, density=False)
    normed, _, _ = np.histogram(data, bins=bins, density=True)

    # Use the middle bin, because it's probably relatively populated
    b = (len(bins) - 1) // 2

    return orig[b] / normed[b]


def _read_ampgen(filename):
    ampgen_branches = (
        ("_1_K~_Px", "_1_K~_Py", "_1_K~_Pz", "_1_K~_E"),
        ("_2_pi#_Px", "_2_pi#_Py", "_2_pi#_Pz", "_2_pi#_E"),
        ("_3_pi#_Px", "_3_pi#_Py", "_3_pi#_Pz", "_3_pi#_E"),
        ("_4_pi~_Px", "_4_pi~_Py", "_4_pi~_Pz", "_4_pi~_E"),
    )
    return 1000 * reweight_utils.read_invariant_masses(
        filename, "DalitzEventList", *ampgen_branches
    )


def _read_mc(filename):
    """
    Throws away events with bkgcat != 0

    """
    tree = "TupleDstToD0pi_D0ToKpipipi/DecayTree"
    mc_branches = (
        ("K_PX", "K_PY", "K_PZ", "K_PE"),
        ("pi1_PX", "pi1_PY", "pi1_PZ", "pi1_PE"),
        ("pi2_PX", "pi2_PY", "pi2_PZ", "pi2_PE"),
        ("pi3_PX", "pi3_PY", "pi3_PZ", "pi3_PE"),
    )
    data = reweight_utils.read_invariant_masses(filename, tree, *mc_branches)

    # bkgcat = reweight_utils.read_branch(filename, tree, "Dstar_BKGCAT")
    # delete_indices = [i for i, x in enumerate(bkgcat) if x != 0]
    # np.delete(data, delete_indices, axis=0)

    return data


def _efficiencies(bins, data, model):
    """
    Find the efficiency projections data/model.

    Bins the data into histograms, scaling them to have equal area.

    Returns a 5-element iterable of efficiencies, and a 5-element iterable of their errors

    """
    dim = len(data[0])
    assert len(model[0]) == dim, "Model and data dimensionality"
    assert len(bins) == dim, "Number of bin arrays"

    scaled_binned_data, scaled_binned_model = [None] * dim, [None] * dim
    ratios, errors = [None] * dim, [None] * dim
    for i in range(dim):
        scaled_binned_data[i], _ = np.histogram(data[:, i], bins=bins[i], density=True)
        scaled_binned_model[i], _ = np.histogram(
            model[:, i], bins=bins[i], density=True
        )

        # Find the ratios of model/data
        ratios[i] = np.divide(scaled_binned_data[i], scaled_binned_model[i])

        # Find the error on this ratio
        # Error = scaled ratio * fractional error in unscaled ratio
        errors[i] = np.multiply(
            ratios[i],
            reweight_utils.fractional_ratio_error(bins[i], data[:, i], model[:, i]),
        )

    return ratios, errors


def _plot_hists(mc, model, name, path):
    plt.figure()
    gs.GridSpec(2, 3)

    # First plot as title
    plt.subplot2grid((2, 3), (0, 0))
    plt.text(0.2, 0.5, name, fontsize=16, bbox={"color": "white"})
    plt.axis("off")

    kw = {"bins": 75, "alpha": 0.4, "density": True}
    for i, ax in enumerate(([0, 1], [0, 2], [1, 0], [1, 1], [1, 2])):
        plt.subplot2grid((2, 3), ax)
        plt.hist(model[:, i], label="Model", **kw)

    plt.savefig(path)


def _flat_efficiency(bins):
    print("Reading flat data")
    flat_data = _read_mc("2018MCflat.root")

    # Remember to perform momentum ordering
    k, pi1, pi2, pi3 = script_util.flat_phsp_points(len(flat_data))
    for i in range(len(k[0])):
        pi1.T[i], pi2.T[i] = reweight_utils.momentum_order(k.T[i], pi1.T[i], pi2.T[i])

    flat_model = reweight_utils.invariant_mass_parametrisation(k, pi1, pi2, pi3)

    _plot_hists(flat_data, flat_model, "flat", "flat.png")
    return _efficiencies(bins, flat_data, flat_model)


def _ampgen_efficiency(bins, mc_filename, ampgen_filename):
    print(f"Reading {ampgen_filename}")
    ampgen_data = _read_ampgen(ampgen_filename)

    print(f"Reading {mc_filename}")
    mc_data = _read_mc(mc_filename)
    _plot_hists(mc_data, ampgen_data, mc_filename, f"{mc_filename}_proj.png")

    return _efficiencies(bins, mc_data, ampgen_data)


def main():
    # Define bins
    n = 100
    bins = [
        np.linspace(600, 1200, n),
        np.linspace(350, 1200, n),
        np.linspace(350, 1200, n),
        np.linspace(900, 1800, n),
        np.linspace(550, 1400, n),
    ]
    bin_centres = [(b[1:] + b[:-1]) / 2 for b in bins]

    # Calculate the flat, RS and WS efficiency projections
    flat = _flat_efficiency(bins)
    # rs = _ampgen_efficiency(bins, "2018MC_RS.root", "ampgen_Dbar_RS.root")
    # ws = _ampgen_efficiency(bins, "2018MC_WS.root", "ampgen_WS.root")

    # Plot them
    plotting.plot_efficiencies(
        bin_centres,
        flat,
        # rs, ws,
        labels=[
            "flat",
            # "rs", "ws
        ],
        path="tmp.png",
    )


if __name__ == "__main__":
    main()
