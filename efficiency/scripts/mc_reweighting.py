import numpy as np
import matplotlib.pyplot as plt
import phasespace
import sys
import os
from sklearn.model_selection import train_test_split

sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
import classification
import reweight_utils
import reweighting
import script_util


def where(file_name: str, tree_name: str, branch_name: str, f):
    """
    Find indices on the specified branch where f(x) = True

    """
    data = reweight_utils.read_branch(file_name, tree_name, branch_name)

    # Might be slow for large arrays. Could use itertools.compress or numpy
    return [i for i, x in enumerate(data) if f(x)]


def flat_k3pi_events(num_events):
    """
    Generate a lot of flat k3pi events.

    Returns a numpy array of shape (num_events, 5) where each 5-element entry is the phase space point in the
    invariant mass parameterisation

    """
    # Charged pion/kaon, neutral D. These might be wrong
    pi_mass = 139.570
    k_mass = 493.677
    d_mass = 1864.84

    # Find the kinematic information
    # Could possibly make this more efficient by doing it in a loop, maybe
    _, particles = phasespace.nbody_decay(
        d_mass, (k_mass, pi_mass, pi_mass, pi_mass)
    ).generate(n_events=num_events)

    # Create labels for the particles to make it easier to code because im tired
    k = particles["p_0"].numpy()
    pi1 = particles["p_1"].numpy()
    pi2 = particles["p_2"].numpy()
    pi3 = particles["p_3"].numpy()

    # Translate this into phase space points:
    #     m(Kpi1), m(pi1pi2), m(pi2pi3), m(kpi1pi2), m(pi1pi2pi3)
    return reweight_utils.invariant_mass_parametrisation(k.T, pi1.T, pi2.T, pi3.T)


def mc_k3pi_events(file_name, tree_name):
    """
    Read LHCb MC events from a file that I've prepared

    """
    mc_events = reweight_utils.read_invariant_masses(
        file_name,
        tree_name,
        ("K_PX", "K_PY", "K_PZ", "K_PE"),
        ("pi1_PX", "pi1_PY", "pi1_PZ", "pi1_PE"),
        ("pi2_PX", "pi2_PY", "pi2_PZ", "pi2_PE"),
        ("pi3_PX", "pi3_PY", "pi3_PZ", "pi3_PE"),
    )

    # Delete all non-signal events
    delete_indices = where(file_name, tree_name, "Dstar_BKGCAT", lambda x: x != 0)
    return np.delete(mc_events, delete_indices, axis=0)


def plot_projections(phsp, mc, weights):
    """
    Plot phsp projections of phsp and un/weighted MC data

    """
    xlabels = (
        r"$M(K\pi_1)$",
        r"$M(\pi_1\pi_2)$",
        r"$M(\pi_2\pi_3)$",
        r"$M(K\pi_1\pi_2)$",
        r"$M(\pi_1\pi_2\pi_3)$",
    )
    shape = (2, 6)
    loc = ((0, 0), (0, 2), (0, 4), (1, 1), (1, 3))
    bins = np.linspace(0, 2000, 125)
    for i in range(5):
        # We want 5 subplots for the figure
        plt.subplot2grid(shape, loc[i], colspan=2)

        # Plot our flat events with no weighting
        plt.hist(phsp[:, i], bins=bins, alpha=0.3, label="phsp", edgecolor="k")

        # Plot LHCb MC events with no weighting
        plt.hist(mc[:, i], bins=bins, alpha=0.3, label="MC", edgecolor="k")

        # Plot the LHCb MC events with the BDT reweighting
        plt.hist(
            mc[:, i],
            weights=weights,
            bins=bins,
            alpha=0.3,
            label="Reweighted MC",
            edgecolor="k",
        )

        # Only plot the legend on the first subplot
        if not i:
            plt.legend(loc="upper left")

        plt.ylabel("Counts")
        plt.xlabel(xlabels[i])

    plt.suptitle("MC Phsp Projections")
    plt.show()


def plot_diffs(phsp, mc, weights):
    """
    Plot phsp projections along side their delta from MC

    """
    # Need to bin data first
    kwargs = {"bins": 100, "range": (0, 2000)}
    xlabels = (
        r"$M(K\pi_1)$",
        r"$M(\pi_1\pi_2)$",
        r"$M(\pi_2\pi_3)$",
        r"$M(K\pi_1\pi_2)$",
        r"$M(\pi_1\pi_2\pi_3)$",
    )
    for i in range(5):
        mc_hist, bins = np.histogram(mc[:, i], **kwargs)
        phsp_hist, _ = np.histogram(phsp[:, i], **kwargs)
        mc_weighted, _ = np.histogram(mc[:, i], weights=weights, **kwargs)

        # Find bin centres and widths
        centres = [(bins[i + 1] + bins[i]) / 2 for i in range(len(bins) - 1)]
        widths = [bins[i + 1] - bins[i] for i in range(len(bins) - 1)]

        script_util.plot_hist_diffs(
            mc_hist,
            mc_weighted,
            phsp_hist,
            centres,
            np.divide(widths, 2),
            xlabels[i],
            ("MC", "Reweighted", "Flat"),
            ("MC-Flat", "Reweighted-Flat"),
            "Phsp Projection",
        )


def roc_curve(mc, flat, weights):
    """
    Show a plot of the ROC curves for a simple binary classifier, both before + after applying weights to the MC data

    """
    curve_before = classification.roc_curve(
        mc, flat, np.ones_like(mc[:, 0]), np.ones_like(flat[:, 0])
    )
    curve_after = classification.roc_curve(mc, flat, weights, np.ones_like(flat[:, 0]))

    plt.plot(curve_before[0], curve_before[1], label="Before")
    plt.plot(curve_after[0], curve_after[1], label="After")
    linspace = np.linspace(0, 1)
    plt.plot(linspace, linspace, label="Indistinguishable", linestyle="--", color="k")

    plt.legend()
    plt.xlabel("FPR")
    plt.ylabel("TPR")
    plt.title("ROC for classification of MC and phsp events")

    plt.show()


def main():
    file_name = "2018MCflat.root"
    tree_name = "TupleDstToD0pi_D0ToKpipipi/DecayTree"

    # Read phsp events from file
    mc_events = mc_k3pi_events(file_name, tree_name)

    # Generate a load of flat d->k3pi events in phase space
    flat_events = flat_k3pi_events(len(mc_events))

    # Split data into training and testing samples
    mc_train, mc_test, flat_train, flat_test = train_test_split(mc_events, flat_events)

    # Train the reweighting BDT
    print("Training BDT...")
    bdt = reweighting.init(flat_train, mc_train, learning_rate=0.15, n_estimators=80)

    # Find the weights from this BDT for our training data
    print("Finding weights...")
    weights = reweighting.predicted_weights(
        bdt, mc_test, expected_num_points=len(mc_test)
    )

    # Plot projections
    plot_projections(flat_test, mc_test, weights)

    # Plot diffs
    plot_diffs(flat_test, mc_test, weights)

    # Plot a ROC curve
    print("Calculating ROC curve...")
    roc_curve(mc_test, flat_test, weights)


if __name__ == "__main__":
    main()
