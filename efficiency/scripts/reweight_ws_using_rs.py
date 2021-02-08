"""
Reweight the WS distribution using weights calculated from the RS MC and model

"""
import os
import sys
import numpy as np
from matplotlib import use as mpl_use
import matplotlib.pyplot as plt
from sklearn.ensemble import GradientBoostingClassifier

mpl_use("Agg")
sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
sys.path.append(os.path.dirname(__file__) + "/../goodness_of_fit/")
import classification
import script_util
import reweight_utils
import goodness_of_fit
import plotting
import reweighting


def _read_ampgen(filename):
    ampgen_branches = (
        ("_1_K~_Px", "_1_K~_Py", "_1_K~_Pz", "_1_K~_E"),
        ("_2_pi#_Px", "_2_pi#_Py", "_2_pi#_Pz", "_2_pi#_E"),
        ("_3_pi#_Px", "_3_pi#_Py", "_3_pi#_Pz", "_3_pi#_E"),
        ("_4_pi~_Px", "_4_pi~_Py", "_4_pi~_Pz", "_4_pi~_E"),
    )
    return reweight_utils.read_invariant_masses(
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

    bkgcat = reweight_utils.read_branch(filename, tree, "Dstar_BKGCAT")
    delete_indices = [i for i, x in enumerate(bkgcat) if x != 0]
    np.delete(data, delete_indices, axis=0)

    return data


def main():
    # read WS MC/AmpGen
    print("Reading WS AmpGen")
    ws_ampgen = _read_ampgen("ampgen_WS.root")
    print("Reading WS MC")
    ws_mc = _read_mc("2018MC_WS.root")

    # Train a classifier to distinguish the WS points
    # Label AmpGen points 0 and MC points 1
    labels = np.concatenate((np.zeros(len(ws_ampgen)), np.ones(len(ws_mc))))
    data = np.concatenate((ws_ampgen, ws_mc), axis=0)
    Classifier = GradientBoostingClassifier(n_estimators=200).fit(data, labels)

    # Find how well the distributions agree
    bins = np.linspace(0, 1)
    ampgen_prob, mc_prob, chisq, p = goodness_of_fit.goodness_of_fit(
        ws_mc, ws_ampgen, Classifier, bins
    )

    # Plot classification
    plotting.classification_plots(
        ws_mc,
        ws_ampgen,
        mc_prob,
        ampgen_prob,
        "hello",
        "WS AmpGen",
        chisq,
        p,
        bins,
        path="WS_class.png",
    )

    # Read RS MC/AmpGen, deleting non-signal events
    print("Reading RS AmpGen")
    rs_ampgen = _read_ampgen("ampgen_Dbar_RS.root")
    print("Reading RS MC")
    rs_mc = _read_mc("2018MC_RS.root")

    # train BDT on RS
    bdt = reweighting.init(rs_ampgen, rs_mc, n_estimators=150)

    # Find weights TODO normalise
    mc_weights = reweighting.predicted_weights(bdt, ws_mc)
    weights = np.concatenate((np.ones(len(ws_ampgen)), mc_weights))

    # Re train the classifier
    Classifier = GradientBoostingClassifier(n_estimators=200).fit(data, labels, weights)
    ampgen_prob, mc_prob, chisq, p = goodness_of_fit.goodness_of_fit(
        ws_mc, ws_ampgen, Classifier, bins, original_weights=mc_weights
    )

    plotting.classification_plots(
        ws_mc,
        ws_ampgen,
        mc_prob,
        ampgen_prob,
        "hello",
        "WS AmpGen",
        chisq,
        p,
        bins,
        mc_weights=mc_weights,
        path="WS_class2.png",
    )


if __name__ == "__main__":
    main()
