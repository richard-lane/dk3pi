"""
An example showing the multidimensional goodness of fit using the 5D final-state phase space of a D -> K3pi decay.

Uses phsp data from LHCb MC, and dynamically generated phsp data

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
import goodness_of_fit
import script_util
import reweight_utils
import reweighting


def _generate_flat_data(num_points):
    """
    Generate flat D->K3pi events, uniformly distributed in phase space

    Returns the momentum-ordered phase space parametrisation

    """
    flat_k, flat_pi1, flat_pi2, flat_pi3 = script_util.flat_phsp_points(num_points)

    # Order the flat points so that pi1 and pi2 are in order of M(kpi)
    for i in range(len(flat_k[0])):
        flat_pi1.T[i], flat_pi2.T[i] = reweight_utils.momentum_order(
            flat_k.T[i], flat_pi1.T[i], flat_pi2.T[i]
        )

    # Find invariant masses
    return reweight_utils.invariant_mass_parametrisation(
        flat_k, flat_pi1, flat_pi2, flat_pi3
    )


def _train_classifier(mc, flat, mc_weights=None, flat_weights=None):
    """
    Create and train a GradientBoostingClassifier

    Returns the trained classifier

    """
    if mc_weights is not None and len(mc_weights) != len(mc):
        raise ValueError(
            f"MC data and weights have {len(mc)} and {len(mc_weights)} entries respectively"
        )

    if flat_weights is not None and len(flat) != len(flat_weights):
        raise ValueError(
            f"Flat distribution and weights have {len(flat)} and {len(flat_weights)} entries respectively"
        )

    # Label training and testing data, and concatenate data into one contiguous thing
    # 0 for MC, 1 for flat
    labels = np.concatenate((np.zeros(len(mc)), np.ones(len(flat))))
    data = np.concatenate((mc, flat), axis=0)

    if mc_weights is None and flat_weights is None:
        weights = None
    elif mc_weights is None:
        weights = np.concatenate((np.ones(len(mc_weights)), flat_weights))
    elif flat_weights is None:
        weights = np.concatenate((mc_weights, np.ones(len(flat))))
    else:
        weights = np.concatenate((mc_weights, flat_weights))

    return GradientBoostingClassifier(n_estimators=200).fit(data, labels, weights)


def _plots(mc, flat, mc_prob, flat_prob, title, chisq, p, prob_bins, mc_weights=None):
    """
    Make plots of the phsp projections of our flat and MC data, and also plot the
    classification probabilities

    """
    plt.figure()
    gs.GridSpec(3, 4)
    plt.suptitle(title)

    # Use the first plot as a title
    plt.subplot2grid((3, 4), (0, 0))
    plt.text(0.2, 0.5, "Phase Space\nProjections", fontsize=16, bbox={"color": "white"})
    plt.axis("off")

    # Large plot to show classification probabilities
    kw = {"alpha": 0.3}
    plt.subplot2grid((3, 4), (0, 2), colspan=2, rowspan=3)
    plt.hist(mc_prob, **kw, label="MC", bins=prob_bins, weights=mc_weights)
    plt.hist(flat_prob, **kw, label="Flat", bins=prob_bins)
    plt.legend()
    plt.yticks([])
    plt.title("Classification probability")
    plt.text(0.1, 0.8, f"$\chi^2$={chisq}\np={p}", transform=plt.gca().transAxes)

    # Small plots to show phsp projections
    for i, ax, label in zip(
        range(5),
        ([0, 1], [1, 0], [1, 1], [2, 0], [2, 1]),
        (
            r"M($K\pi_1$)",
            r"M($\pi_1\pi_2$)",
            r"M($\pi_2\pi_3$)",
            r"M($K\pi_1\pi_2$)",
            r"M($\pi_1\pi_2\pi_3$)",
            r" ",
        ),
    ):
        plt.subplot2grid((3, 4), ax)
        plt.hist(mc[:, i], **kw, label="MC", weights=mc_weights, bins=75)
        plt.hist(flat[:, i], **kw, label="Flat", bins=75)
        plt.yticks([])
        plt.title(label)

        # Only need one legend
        if not i:
            plt.legend()

        # Bottom two need axis labels
        if i in {3, 4}:
            plt.xlabel("MeV")

    plt.show()


def main():
    # Read in LHCb MC data
    filename = "2018MCflat.root"
    tree_name = "TupleDstToD0pi_D0ToKpipipi/DecayTree"
    mc = reweight_utils.read_invariant_masses(
        filename,
        tree_name,
        ("K_PX", "K_PY", "K_PZ", "K_PE"),
        ("pi1_PX", "pi1_PY", "pi1_PZ", "pi1_PE"),
        ("pi2_PX", "pi2_PY", "pi2_PZ", "pi2_PE"),
        ("pi3_PX", "pi3_PY", "pi3_PZ", "pi3_PE"),
    )

    def cut(point):
        return (
            not 630 < point[0] < 1180
            or not 320 < point[1] < 1150
            or not 320 < point[2] < 1150
            or not 870 < point[3] < 1700
            or not 570 < point[4] < 1350
        )

    # Throw away MC where the bkgcat isnt 0 or where the ipchi2 is bad
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

    # Generate flat data
    print(f"{len(mc)} MC events")
    flat = _generate_flat_data(len(mc))

    # Throw away some flat points
    to_delete = np.array([], dtype=int)
    for i, point in enumerate(flat):
        if cut(point):
            to_delete = np.append(to_delete, i)
    flat = np.delete(flat, to_delete, axis=0)

    # Split data into test + train
    # This training data is used to both train the reweighting BDT and the binary classifier
    mc_train, mc_test = train_test_split(mc)
    flat_train, flat_test = train_test_split(flat)

    # Generate uniform weights to apply to the MC data so that our histograms have the same areas
    _, num_mc_test = len(mc_train), len(mc_test)
    _, num_flat_test = len(flat_train), len(flat_test)
    weights = np.ones(num_mc_test) * (num_flat_test / num_mc_test)

    # Train the classifier on unweighted data
    print("Training Classifier")
    Classifier = _train_classifier(mc_train, flat_train)

    # Find their goodness of fit and plot
    bins = np.concatenate(
        ([0.0, 0.35, 0.38], np.linspace(0.385, 0.61, num=50), [0.63, 1.0])
    )
    mc_prob, flat_prob, chisq, p = goodness_of_fit.goodness_of_fit(
        mc_test, flat_test, Classifier, bins, original_weights=weights
    )
    _plots(
        mc_test,
        flat_test,
        mc_prob,
        flat_prob,
        "Unweighted Phsp MC",
        chisq,
        p,
        bins,
        mc_weights=weights,
    )

    # Reweight the MC data to flat
    print("Training BDT")
    bdt = reweighting.init(flat_train, mc_train)
    weights = reweighting.predicted_weights(bdt, mc_test)

    # Scale the weights so that our histograms have the same areas
    weights *= num_flat_test / (num_mc_test * np.mean(weights))

    print("Training Classifier")
    Classifier = _train_classifier(mc_train, flat_train)

    mc_prob, flat_prob, chisq, p = goodness_of_fit.goodness_of_fit(
        mc_test, flat_test, Classifier, bins, weights
    )
    _plots(
        mc_test,
        flat_test,
        mc_prob,
        flat_prob,
        "weighted Phsp MC",
        chisq,
        p,
        bins,
        mc_weights=weights,
    )


if __name__ == "__main__":
    main()
