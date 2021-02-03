import numpy as np
import phasespace
import sys
import os
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from multiprocessing import Process
import matplotlib.gridspec as gs
from tqdm import tqdm

# Don't attempt to create figure windows, since this will probably get run on lxplus
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
sys.path.append(os.path.dirname(__file__) + "/../goodness_of_fit/")
import classification
import reweight_utils
import reweighting
import script_util
import goodness_of_fit


def _train_classifier(mc, model, mc_weights=None, model_weights=None):
    """
    Create and train a GradientBoostingClassifier

    Returns the trained classifier

    """
    if mc_weights is not None and len(mc_weights) != len(mc):
        raise ValueError(
            f"MC data and weights have {len(mc)} and {len(mc_weights)} entries respectively"
        )

    if model_weights is not None and len(model_weights) != len(model):
        raise ValueError(
            f"Model distribution and weights have {len(model)} and {len(model_weights)} entries respectively"
        )

    # Label training and testing data, and concatenate data into one contiguous thing
    # 0 for MC, 1 for model
    labels = np.concatenate((np.zeros(len(mc)), np.ones(len(model))))
    data = np.concatenate((mc, model), axis=0)

    if mc_weights is None and model_weights is None:
        weights = None
    elif mc_weights is None:
        weights = np.concatenate((np.ones(len(mc_weights)), model_weights))
    elif flat_weights is None:
        weights = np.concatenate((mc_weights, np.ones(len(model_weights))))
    else:
        weights = np.concatenate((mc_weights, model_weights))

    return GradientBoostingClassifier(n_estimators=200).fit(data, labels, weights)


def _plots(
    mc,
    model,
    mc_prob,
    model_prob,
    title,
    chisq,
    p,
    prob_bins,
    path,
    label,
    mc_weights=None,
):
    """
    Make plots of the phsp projections of our model and MC data, and also plot the
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
    plt.hist(model_prob, **kw, label=label, bins=prob_bins)
    plt.legend()
    plt.yticks([])
    plt.title("Classification probability")
    plt.text(0.1, 0.9, f"$\chi^2$={chisq}\np={p}")

    # Small plots to show phsp projections
    for i, ax, ax_label in zip(
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
        plt.hist(model[:, i], **kw, label=label, bins=75)
        plt.yticks([])
        plt.title(ax_label)

        # Only need one legend
        if not i:
            plt.legend()

        # Bottom two need axis labels
        if i in {3, 4}:
            plt.xlabel("MeV")

    plt.savefig(path)


def where(file_name: str, tree_name: str, branch_name: str, f):
    """
    Find indices on the specified branch where f(x) = True

    Not very fast and possibly something already exists in numpy

    """
    data = reweight_utils.read_branch(file_name, tree_name, branch_name)

    # Might be slow for large arrays. Could use itertools.compress or numpy
    return [i for i, x in enumerate(data) if f(x)]


def train_bdt(target, origin):
    """
    Helper fcn

    """
    return reweighting.init(target, origin, n_estimators=100, learning_rate=0.05)


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


def plot_diffs(target, mc, weights, path, label):
    """
    Save plots of phsp projections alongside their delta from MC

    Normalises the target disribution to have the same number of counts as mc

    Saves to "diff_{path}_1.png", etc.

    """
    # Find what normalisation we need to apply to the target histograms
    normalisation = len(mc) / len(target)

    kwargs = {"bins": 100, "range": (0, 2000)}
    xlabels = (
        r"$M(K\pi_1)$",
        r"$M(\pi_1\pi_2)$",
        r"$M(\pi_2\pi_3)$",
        r"$M(K\pi_1\pi_2)$",
        r"$M(\pi_1\pi_2\pi_3)$",
    )
    for i in range(5):
        # Need to bin data
        mc_hist, bins = np.histogram(mc[:, i], **kwargs)
        target_hist, _ = np.histogram(target[:, i], **kwargs)

        # Normalise the target hist to the right number of events
        target_hist = target_hist * normalisation

        # Normalise the reweighted hist to the right number of events
        n_reweighted = len(mc) * np.mean(weights)
        n_target = np.sum(target_hist)
        weights *= n_target / n_reweighted

        mc_weighted, _ = np.histogram(mc[:, i], weights=weights, **kwargs)

        # Find bin centres and widths
        centres = [(bins[i + 1] + bins[i]) / 2 for i in range(len(bins) - 1)]
        widths = [bins[i + 1] - bins[i] for i in range(len(bins) - 1)]

        script_util.plot_hist_diffs(
            mc_hist,
            mc_weighted,
            target_hist,
            centres,
            np.divide(widths, 2),
            xlabels[i],
            ("MC", "Reweighted", label),
            (f"MC-{label}", f"Reweighted-{label}"),
            f"{label} Phsp Projection",
        )
        plt.savefig(f"diff_{path}_{i}.png", dpi=600)


def roc_curve(mc, flat, weights, label):
    """
    Create a plot of the ROC curves for a simple binary classifier, both before + after applying weights to the MC data

    """
    plt.clf()
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
    plt.title(f"ROC for classification of MC and {label} events")
    plt.savefig(f"{label}.png")


def flat_study():
    """
    Plot graphs, do reweighting etc. with flat phsp

    """
    file_name = "2018MCflat.root"
    tree_name = "TupleDstToD0pi_D0ToKpipipi/DecayTree"

    # Read phsp events from file
    print("Reading phsp MC file")
    phsp_mc_events = mc_k3pi_events(file_name, tree_name)

    # Generate a load of flat d->k3pi events in phase space
    print("Generating phsp events")
    phsp_flat_events = flat_k3pi_events(len(phsp_mc_events))

    # Delete the events with m(kPi1) > 1300
    print("Deleting events...")
    max_mkpi = 1200
    phsp_flat_events = np.delete(
        phsp_flat_events, phsp_flat_events[:, 0] > max_mkpi, axis=0
    )
    phsp_mc_events = np.delete(phsp_mc_events, phsp_mc_events[:, 0] > max_mkpi, axis=0)

    # Split data into training and testing samples
    (
        mc_train,
        mc_test,
    ) = train_test_split(phsp_mc_events)
    flat_train, flat_test = train_test_split(phsp_flat_events)

    # Train the reweighting BDT
    print("Training BDT...")
    bdt = train_bdt(flat_train, mc_train)

    # Find the weights from this BDT for our test data
    print("Finding weights...")
    weights = reweighting.predicted_weights(
        bdt, mc_test, expected_num_points=len(mc_test)
    )

    # Plot diffs
    plot_diffs(flat_test, mc_test, weights, "Flat", "Phsp MC")
    plot_diffs(
        flat_train,
        mc_train,
        reweighting.predicted_weights(bdt, mc_train, expected_num_points=len(mc_test)),
        "Flat_train",
        "phsp",
    )

    # Plot a ROC curve
    print("Calculating ROC curve...")
    roc_curve(mc_test, flat_test, weights, "Flat")


def ampgen_study(mc_file_name, tree_name, ampgen_file_name, label):
    """
    Helper fcn for doing the RS/WS study with AmpGen data

    """
    # Ugly way of storing params
    branches = (
        ("_1_K~_Px", "_1_K~_Py", "_1_K~_Pz", "_1_K~_E"),
        ("_2_pi#_Px", "_2_pi#_Py", "_2_pi#_Pz", "_2_pi#_E"),
        ("_3_pi#_Px", "_3_pi#_Py", "_3_pi#_Pz", "_3_pi#_E"),
        ("_4_pi~_Px", "_4_pi~_Py", "_4_pi~_Pz", "_4_pi~_E"),
    )

    # Read MC events from file
    print(f"{label}: Reading MC...")
    mc_events = mc_k3pi_events(mc_file_name, tree_name)

    # Read AmpGen events from file
    print(f"{label}: Reading AmpGen...")
    ampgen_events = reweight_utils.read_invariant_masses(
        ampgen_file_name,
        "DalitzEventList",
        branches[0],
        branches[1],
        branches[2],
        branches[3],
    )

    # Convert AmpGen events to MeV
    ampgen_events *= 1000

    # Throw away events with m(kpi1) > 900
    # print("Deleting events...")
    # ampgen_events = np.delete(ampgen_events, ampgen_events[:, 0] > 900, axis=0)
    # mc_events = np.delete(mc_events, mc_events[:, 0] > 900, axis=0)

    # Split data into train/test
    kwargs = {"test_size": 0.5}
    mc_train, mc_test = train_test_split(mc_events, **kwargs)
    ampgen_train, ampgen_test = train_test_split(ampgen_events, **kwargs)

    # Create a classifier to distinguish MC and AmpGen data
    Classifier = _train_classifier(mc_train, ampgen_train)
    # Weights to make the histograms the same area
    norm_weights = np.ones(len(mc_train)) * (len(ampgen_train) / len(mc_train))
    bins = np.concatenate(([0.0], np.linspace(0.2, 0.8, num=10), [1.0]))
    mc_prob, ampgen_prob, chisq, p = goodness_of_fit.goodness_of_fit(
        mc_test, ampgen_test, Classifier, bins, norm_weights
    )

    # Plot something
    _plots(
        mc_train,
        ampgen_train,
        mc_prob,
        ampgen_prob,
        f"Unweighted {label} prob",
        chisq,
        p,
        bins,
        f"unweighted_{label}.png",
        label,
        mc_weights=norm_weights,
    )

    # Train the reweighting BDT
    print(f"{label}: Training bdt...")
    bdt = train_bdt(ampgen_train, mc_train)

    # Find the weights for the testing data, scaling them
    print(f"{label}: Finding weights...")
    weights = reweighting.predicted_weights(bdt, mc_test)
    weights *= len(ampgen_test) / (len(mc_test) * np.mean(weights))

    # Plot something again
    _plots(
        mc_train,
        ampgen_train,
        mc_prob,
        ampgen_prob,
        f"Unweighted {label} prob",
        chisq,
        p,
        bins,
        f"weighted_{label}.png",
        label,
        mc_weights=norm_weights,
    )

    # Plot diffs
    # plot_diffs(ampgen_test, mc_test, weights, f"{label}_test", f"{label} AmpGen")

    # Plot diffs of training data. Just to check
    # plot_diffs(
    #     ampgen_train,
    #     mc_train,
    #     reweighting.predicted_weights(bdt, mc_train),
    #     f"{label}_train",
    #     f"{label} AmpGen",
    # )

    # Plot a ROC curve
    # print(f"{label}: Calculating ROC curve...")
    # roc_curve(mc_test, ampgen_test, weights, label)


def ws_study():
    """
    Plot graphs etc. for RS data

    """
    ampgen_study(
        "2018MC_WS.root",
        "TupleDstToD0pi_D0ToKpipipi/DecayTree",
        "ampgen_WS.root",
        "WS",
    )


def rs_study():
    """
    Plot graphs etc. for RS data

    """
    ampgen_study(
        "2018MC_RS.root",
        "TupleDstToD0pi_D0ToKpipipi/DecayTree",
        "ampgen_Dbar_RS.root",
        "RS",
    )


def mc_studies():
    ws = Process(target=ws_study)
    rs = Process(target=rs_study)

    ws.start()
    rs.start()

    ws.join()
    rs.join()


def main():
    # flat_study()
    mc_studies()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass
