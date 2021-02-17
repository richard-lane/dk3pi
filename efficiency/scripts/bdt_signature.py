"""
Train two reweighting BDTs on two identical datasets- one BDT takes the decay's signature into account, the other doesn't.

Then train two binary classifiers to identify samples of these events, and compare

"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split

from tqdm import tqdm

sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), "scripts"))
)
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__))))
)
sys.path.append(
    os.path.abspath(
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "goodness_of_fit")
    )
)
import goodness_of_fit
import script_util
import reweight_utils
import reweighting
import plotting


def _generate_flat_data_with_sig(num_points):
    """
    Generate flat D->K3pi events, uniformly distributed in phase space

    Returns the momentum-ordered phase space parametrisation, and the signature of each event

    """
    flat_k, flat_pi1, flat_pi2, flat_pi3 = script_util.flat_phsp_points(num_points)

    # Order the flat points so that pi1 and pi2 are in order of M(kpi)
    for i in range(len(flat_k[0])):
        flat_pi1.T[i], flat_pi2.T[i] = reweight_utils.momentum_order(
            flat_k.T[i], flat_pi1.T[i], flat_pi2.T[i]
        )

    # Find signatures
    sigs = reweight_utils.signatures(flat_k, flat_pi1, flat_pi2, flat_pi3)

    # Find invariant masses
    return (
        reweight_utils.invariant_mass_parametrisation(
            flat_k, flat_pi1, flat_pi2, flat_pi3
        ),
        sigs,
    )


def _read_mc_with_sig(filename, tree_name):
    k_branches, pi1_branches, pi2_branches, pi3_branches = (
        ("K_PX", "K_PY", "K_PZ", "K_PE"),
        ("pi1_PX", "pi1_PY", "pi1_PZ", "pi1_PE"),
        ("pi2_PX", "pi2_PY", "pi2_PZ", "pi2_PE"),
        ("pi3_PX", "pi3_PY", "pi3_PZ", "pi3_PE"),
    )
    k, pi1, pi2, pi3 = reweight_utils.read_kinematic_data(
        filename, tree_name, k_branches, pi1_branches, pi2_branches, pi3_branches
    )

    for i in range(len(k[0])):
        # Assign i'th pi1 and pi2 params to the right things
        pi1.T[i], pi2.T[i] = reweight_utils.momentum_order(k.T[i], pi1.T[i], pi2.T[i])

    mc = reweight_utils.invariant_mass_parametrisation(k, pi1, pi2, pi3)
    sig = reweight_utils.signatures(k, pi1, pi2, pi3)
    return mc, sig


def _analyse(data, model, bins, title):
    """
    Split data/model into test+train, train a BDT to reweight, train a binary classifier to separate them, plot the classification probabilities

    """
    # Split into test/train
    print("splitting")
    model_train, model_test, = train_test_split(model, train_size=0.5)
    data_train, data_test, = train_test_split(data, train_size=0.5)

    # Train a BDT, find weights
    print("training bdt")
    bdt = reweighting.init(model_train, data_train, n_estimators=200)
    print("finding weights")
    weights = reweighting.predicted_weights(bdt, data_test)

    # Train a classifier
    data_train_all = np.concatenate((data_train, model_train))
    labels = np.concatenate((np.zeros(len(data_train)), np.ones(len(model_train))))
    print("classifiying")
    Classifier = GradientBoostingClassifier(n_estimators=200).fit(
        data_train_all, labels
    )

    # Use the classifier to classify the reweighted events
    print("finding goodness of fit")
    data_prob, model_prob, chisq, p = goodness_of_fit.goodness_of_fit(
        data_test, model_test, Classifier, bins, original_weights=weights
    )
    plotting.classification_plots(
        data_test,
        model_test,
        data_prob,
        model_prob,
        title,
        "Phsp",
        chisq,
        p,
        bins,
        weights,
    )


def main():
    # Read in LHCb MC data
    filename = "2018MCflat.root"
    tree_name = "TupleDstToD0pi_D0ToKpipipi/DecayTree"
    mc, mc_sig = _read_mc_with_sig(filename, tree_name)

    def cut(point):
        return (
            not 630 < point[0] < 1180
            or not 320 < point[1] < 1150
            or not 320 < point[2] < 1150
            or not 870 < point[3] < 1700
            or not 570 < point[4] < 1350
        )

    # Generate flat data
    flat, flat_sig = _generate_flat_data_with_sig(len(mc))

    # Throw away MC where the bkgcat isnt 0 or where the ipchi2 is bad or where the cuts don't pass
    dst_ipchi2 = reweight_utils.read_branch(filename, tree_name, "Dstar_IPCHI2_OWNPV")
    d_ipchi2 = reweight_utils.read_branch(filename, tree_name, "D_IPCHI2_OWNPV")
    bkgcat = reweight_utils.read_branch(filename, tree_name, "Dstar_BKGCAT")
    indices = [
        i
        for i, (bkg, dst_ip, d_ip, event) in enumerate(
            zip(bkgcat, dst_ipchi2, d_ipchi2, mc)
        )
        if bkg != 0 or not 0 < dst_ip < 9 or not 0 < d_ip < 9 or cut(event)
    ]
    mc = np.delete(mc, indices, axis=0)
    mc_sig = np.delete(mc_sig, indices)

    # Throw away flat data where the cuts don't pass
    to_delete = np.array([], dtype=int)
    for i, point in enumerate(flat):
        if cut(point):
            to_delete = np.append(to_delete, i)
    flat = np.delete(flat, to_delete, axis=0)
    flat_sig = np.delete(flat_sig, to_delete)

    # Throw away a random set of flat data until we have the same number of flat and MC points
    gen = np.random.default_rng()
    to_delete = gen.choice(len(flat), size=(len(flat) - len(mc)), replace=False)
    flat = np.delete(flat, to_delete, axis=0)
    flat_sig = np.delete(flat_sig, to_delete)

    bins = np.concatenate(
        ([0.0, 0.22, 0.26], np.linspace(0.3, 0.61, num=50), [0.63, 1.0])
    )
    _analyse(mc, flat, bins, "no sig")

    # Add the signature to the first phsp point as a multiplier
    mc[:, 0] = np.multiply(mc[:, 0], mc_sig)
    flat[:, 0] = np.multiply(flat[:, 0], flat_sig)
    _analyse(mc, flat, bins, "sig")


if __name__ == "__main__":
    main()
