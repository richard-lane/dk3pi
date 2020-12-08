"""
BDT reweighting in pure python

"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import warnings

import skopt
import joblib
from sklearn.metrics import roc_curve

# This is really horrible but i can't think of a better way of doing it
# Ideally i'd like to set a global python include path via the CMake build system...
sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
import reweighting
import reweight_utils
import classification
import visualisations


def chi_squared(counts_source, counts_target):
    """
    Implementation of chisq test for two distributions that ignores zeroes

    Pass in counts in each bin

    ChiSq defined as Sum((x_i - y_i)^2/(x_i + y_i)), because I do what I want

    This should be unit tested but it isnt. But it works

    """
    assert len(counts_source) == len(counts_target)

    # Want to catch warnings as if they were errors
    warnings.filterwarnings("error")

    try:
        chi_sq = 0.0
        for source, target in zip(counts_source, counts_target):
            try:
                x = (target - source) ** 2 / (target + source)
                chi_sq += x

            except RuntimeWarning:
                # Misses the case where target = -source. But that won't happen...
                pass

    finally:
        # Put the warning status back
        warnings.filterwarnings("default")

    return chi_sq


def combined_chi_squared(
    source_points, source_weights, target_points, target_weights, bins
):
    """
    Take collections of multidimensional points + (scalar) weights and a binning to apply to each dimension

    Then finds the chi squared of each histogram projection using this binning and adds them to return a combined chi squared value

    """
    dimensionality = len(source_points[0])
    assert dimensionality == len(target_points[0])

    chi_sq = 0.0

    for d in range(dimensionality):

        # Bin the points from this dimension into a histogram
        source_counts, edges = np.histogram(
            source_points[:, d], bins=bins, weights=source_weights
        )
        target_counts, _ = np.histogram(
            target_points[:, d], bins=bins, weights=target_weights
        )

        assert np.all(
            np.abs(edges - bins) < 1e-3
        )  # Just in case it isn't. Luke told me to watch out

        # Find chi squared between these histograms
        chi_sq += chi_squared(source_counts, target_counts)

    return chi_sq


def bin_data(
    source, target, reweighted, source_weights, target_weights, reweighted_weights, bins
):
    """
    Returns source counts, target counts, reweighted counts

    """
    source_counts, _ = np.histogram(source, bins=bins, weights=source_weights)
    target_counts, _ = np.histogram(target, bins=bins, weights=target_weights)
    reweighted_counts, _ = np.histogram(
        reweighted, bins=bins, weights=reweighted_weights
    )

    return source_counts, target_counts, reweighted_counts


def rescale(counts, errors):
    """
    Rescale a histogram and its errors to a total area of 1

    """
    assert len(errors) == len(counts)

    integral = np.sum(counts)
    counts /= integral
    errors /= integral


def save_plot(
    bin_centres,
    bin_widths,
    prompt_counts,
    prompt_err,
    reweighted_counts,
    reweighted_err,
    sl_counts,
    sl_err,
    path,
    x_label,
    title,
    plot_errs=True,
):
    """
    Save a plot of our histograms

    """
    fig, ax = plt.subplots(2, sharex=True, gridspec_kw={"height_ratios": [2.5, 1]})
    line_width = 0.5
    markersize = 0.0
    alpha = 0.5
    marker = "_"
    for count, err, colour, label in zip(
        (prompt_counts, reweighted_counts, sl_counts),
        (prompt_err, reweighted_err, sl_err),
        ("red", "blue", "green"),
        ("Prompt", "Reweighted", "SL"),
    ):
        ax[0].errorbar(
            bin_centres,
            count,
            yerr=err if plot_errs else None,
            xerr=bin_widths if plot_errs else None,
            label=label,
            alpha=alpha,
            color=colour,
            linewidth=line_width,
            marker=marker,
            markersize=markersize,
            fmt=" ",
        )

    # Find difference between Prompt + SL and difference between Reweighted + SL
    prompt_minus_sl = prompt_counts - sl_counts
    reweighted_minus_sl = reweighted_counts - sl_counts

    # Find the associated errors
    prompt_minus_sl_err = np.sqrt(prompt_err ** 2 + sl_err ** 2)
    reweighted_minus_sl_err = np.sqrt(reweighted_err ** 2 + sl_err ** 2)

    # Plot
    for diff, err, colour, label in zip(
        (prompt_minus_sl, reweighted_minus_sl),
        (prompt_minus_sl_err, reweighted_minus_sl_err),
        ("red", "blue"),
        ("Prompt-SL", "Reweighted-SL"),
    ):
        ax[1].errorbar(
            bin_centres,
            diff,
            yerr=err if plot_errs else None,
            xerr=bin_widths if plot_errs else None,
            label=label,
            color=colour,
            linewidth=line_width,
            marker=marker,
            markersize=markersize,
            fmt=" ",
        )

    ax[0].legend()
    ax[1].legend()
    ax[0].set_ylabel("Counts (normalised)")
    ax[1].set_ylabel(r"$\Delta$Counts")
    plt.xlabel(x_label)
    fig.subplots_adjust(hspace=0)
    fig.suptitle(title)

    plt.savefig(path, dpi=600, bbox_inches="tight")
    plt.clf()


def read_data():
    """
    Returns phsp parametrised prompt_points, prompt_weights, sl_points, sl_weights

    """
    # Find phsp points for prompt + SL datasets
    print("Reading data...")
    prompt_points = reweight_utils.inv_mass_parametrisation(
        "cut_wg_rs_prompt.root",
        "DecayTree",
        ("D0_P0_PX", "D0_P0_PY", "D0_P0_PZ", "D0_P0_PE"),
        ("D0_P1_PX", "D0_P1_PY", "D0_P1_PZ", "D0_P1_PE"),
        ("D0_P2_PX", "D0_P2_PY", "D0_P2_PZ", "D0_P2_PE"),
        ("D0_P3_PX", "D0_P3_PY", "D0_P3_PZ", "D0_P3_PE"),
    )
    sl_points = reweight_utils.inv_mass_parametrisation(
        "cut_wg_rs_sl.root",
        "DecayTree",
        ("D0_P0_PX", "D0_P0_PY", "D0_P0_PZ", "D0_P0_PE"),
        ("D0_P1_PX", "D0_P1_PY", "D0_P1_PZ", "D0_P1_PE"),
        ("D0_P2_PX", "D0_P2_PY", "D0_P2_PZ", "D0_P2_PE"),
        ("D0_P3_PX", "D0_P3_PY", "D0_P3_PZ", "D0_P3_PE"),
    )

    # Find weights for prompt + SL datasets
    print("Reading weights...")
    prompt_weights = reweight_utils.read_branch(
        "rs_weights.root", "DecayTree", "numSignalEvents_sw"
    )
    sl_weights = reweight_utils.read_branch(
        "sl_weights.root", "DecayTree", "numSignalEvents_sw"
    )

    # Remove some points to make the distributions look more different
    indices_to_delete = []
    for i in range(len(sl_points)):
        if sl_points[i][1] < 900 * np.random.random():
            indices_to_delete.append(i)
    sl_points = np.delete(sl_points, indices_to_delete, axis=0)
    sl_weights = np.delete(sl_weights, indices_to_delete)

    return prompt_points, prompt_weights, sl_points, sl_weights


def make_plots(
    test_prompt_data,
    test_sl_data,
    test_prompt_weights,
    test_sl_weights,
    efficiency_weights,
    bins,
):
    """
    Plot phase space projections, save files to 0.png, 1.png, 2.png etc.

    """
    units = (
        r"M($K\pi_1$) /MeV",
        r"M($\pi_1\pi_2$) /MeV",
        r"M($\pi_2\pi_3$) /MeV",
        r"M($K\pi_1\pi_2$) /MeV",
        r"M($\pi_1\pi_2\pi_3$) /MeV",
    )
    titles = [
        r"Projection: M($K\pi_1$)",
        r"Projection: M($\pi_1\pi_2$)",
        r"Projection: M($\pi_1\pi_2\pi_3$)",
        r"Projection: M($K\pi_1\pi_2$)",
        r"Projection: M($\pi_1\pi_2\pi_3$)",
    ]
    # Compare the reweighted test prompt + test SL data by plotting some histograms
    for i in range(len(test_prompt_data[0])):

        # Find the i'th histogram projection
        prompt, sl, reweighted = bin_data(
            test_prompt_data[:, i],
            test_sl_data[:, i],
            test_prompt_data[:, i],
            test_prompt_weights,
            test_sl_weights,
            efficiency_weights,
            bins,
        )

        # Find errors
        prompt_err = np.sqrt(prompt)
        reweighted_err = np.sqrt(reweighted)
        sl_err = np.sqrt(sl)

        # Rescale histograms and errors
        rescale(prompt, prompt_err)
        rescale(sl, sl_err)
        rescale(reweighted, reweighted_err)

        # Find bin centres
        centres = np.mean(np.vstack([bins[0:-1], bins[1:]]), axis=0)

        # Find bin widths
        widths = [0.5 * (j - i) for i, j in zip(bins[:-1], bins[1:])]

        # Plot
        save_plot(
            centres,
            widths,
            prompt,
            prompt_err,
            reweighted,
            reweighted_err,
            sl,
            sl_err,
            f"{i}.png",
            units[i],
            titles[i],
            True,
        )
        print(f"Created {i}")


def plot_projections():
    """
    Plot phsp projections with default BDT

    """
    # Read data from files + perform phsp parametrisation
    prompt_points, prompt_weights, sl_points, sl_weights = read_data()

    # Split data into training + test data
    print("Splitting data...")
    training_prompt_data, test_prompt_data = np.array_split(prompt_points, 2)
    training_sl_data, test_sl_data = np.array_split(sl_points, 2)
    training_prompt_weights, test_prompt_weights = np.array_split(prompt_weights, 2)
    training_sl_weights, test_sl_weights = np.array_split(sl_weights, 2)

    # Train the BDT on training data
    print("Training BDT...")
    bdt = reweighting.init(
        training_sl_data,
        training_prompt_data,
        training_sl_weights,
        training_prompt_weights,
    )

    # Reweight the test prompt data
    print("Reweighting...")
    efficiency_weights = reweighting.predicted_weights(
        bdt, test_prompt_data, test_prompt_weights
    )

    # Plot phsp projections
    bins = np.linspace(200, 1800, 250)
    make_plots(
        test_prompt_data,
        test_sl_data,
        test_prompt_weights,
        test_sl_weights,
        efficiency_weights,
        bins,
    )


def optimise(n_calls):
    """
    Find the optimal BDT parameters for a given dataset...

    Try n_calls BDT configurations; starts with n_calls//5 random points

    """
    # Read data from files + perform phsp parametrisation
    prompt_points, prompt_weights, sl_points, sl_weights = read_data()

    # Split data into training + test data
    print("Splitting data...")
    training_prompt_data, test_prompt_data = np.array_split(prompt_points, 2)
    training_sl_data, test_sl_data = np.array_split(sl_points, 2)
    training_prompt_weights, test_prompt_weights = np.array_split(prompt_weights, 2)
    training_sl_weights, test_sl_weights = np.array_split(sl_weights, 2)

    bins = np.linspace(200, 1800, 250)

    def objective_fcn(args):
        """
        Train BDT, reweight + find chi squared

        """
        (
            n_estimators,
            learning_rate,
            max_depth,
            min_samples_leaf,
            loss_regularization,
        ) = args

        bdt = reweighting.init(
            training_sl_data,
            training_prompt_data,
            training_sl_weights,
            training_prompt_weights,
            n_estimators=n_estimators,
            learning_rate=learning_rate,
            max_depth=max_depth,
            min_samples_leaf=min_samples_leaf,
            loss_regularization=loss_regularization,
        )

        efficiency_weights = reweighting.predicted_weights(
            bdt, test_prompt_data, test_prompt_weights
        )

        return combined_chi_squared(
            test_prompt_data, efficiency_weights, test_sl_data, test_sl_weights, bins
        )

    dimensions = [
        skopt.space.Integer(20, 250),  # Num trees
        skopt.space.Real(0.01, 0.8, prior="log-uniform"),  # Learning Rate
        skopt.space.Integer(3, 10),  # Tree depth
        skopt.space.Integer(20, 750),  # min samples leaf
        skopt.space.Real(1.5, 50.0),  # loss reg
    ]

    result = skopt.gp_minimize(
        objective_fcn,
        dimensions,
        n_calls=n_calls,
        n_random_starts=n_calls // 5,
        verbose=True,
    )
    print(result.x)

    # Save a bar chart of chi squared values
    fig, ax = plt.subplots()
    x = [_ for _ in range(len(result.func_vals))]
    y = result.func_vals
    ax.bar(x, y)
    ax.set_yscale("log")
    ax.set_xlabel("Trial")
    ax.set_ylabel("ChiSq")
    plt.savefig("chiSqs.png")

    # Save a text file of Chi squareds and params
    with open("results.txt", "w") as f:
        for chisq, params in zip(y, result.x_iters):
            f.write(f"{chisq}\t{params}\n")
        f.write("\n")


def roc_score(prompt_points, sl_points, prompt_weights, sl_weights):
    """
    Compute ROC:AUC score

    """
    # Split data up
    points_train, points_test, labels_train, labels_test, weights_train, weights_test = classification.split_data_for_classification(
        prompt_points, sl_points, prompt_weights, sl_weights
    )

    # Train classifier
    classifier = classification.train_classifier(
        points_train, labels_train, weights_train
    )

    # Compute score
    return classification.classification_score(
        classifier, points_test, labels_test, weights_test
    )


def roc_score_test():
    """
    Compare ROC AUC scores for unweighted + reweighted distributions

    """
    # Read in data
    print("Reading data...")
    prompt_points, prompt_weights, sl_points, sl_weights = read_data()

    # Find the appropriate BDT weights
    training_prompt_data, test_prompt_data = np.array_split(prompt_points, 2)
    training_sl_data, test_sl_data = np.array_split(sl_points, 2)
    training_prompt_weights, test_prompt_weights = np.array_split(prompt_weights, 2)
    training_sl_weights, test_sl_weights = np.array_split(sl_weights, 2)
    print("Training BDT...")
    bdt = reweighting.init(
        training_sl_data,
        training_prompt_data,
        training_sl_weights,
        training_prompt_weights,
    )
    efficiency_weights = reweighting.predicted_weights(
        bdt, test_prompt_data, test_prompt_weights
    )

    linspace = np.linspace(0, 1)
    curve_before = classification.roc_curve(
        test_prompt_data, test_sl_data, test_prompt_weights, test_sl_weights
    )
    curve_after = classification.roc_curve(
        test_prompt_data, test_sl_data, efficiency_weights, test_sl_weights
    )
    plt.plot(curve_before[0], curve_before[1], label="Before Reweighting")
    plt.plot(curve_after[0], curve_after[1], label="After Reweighting")
    plt.plot(linspace, linspace, label="Indistinguishable", linestyle="--", color="k")
    plt.legend()
    plt.xlabel("FPR")
    plt.ylabel("TPR")
    plt.title("ROC Curves for Binary Classification of Prompt/Semileptonic Phsp Data")
    plt.show()


def plot_slices():
    # Read data
    prompt_points, prompt_weights, sl_points, sl_weights = read_data()
    training_prompt_data, test_prompt_data = np.array_split(prompt_points, 2)
    training_sl_data, test_sl_data = np.array_split(sl_points, 2)
    training_prompt_weights, test_prompt_weights = np.array_split(prompt_weights, 2)
    training_sl_weights, test_sl_weights = np.array_split(sl_weights, 2)

    # Reweight
    print("Training BDT...")
    bdt = reweighting.init(
        training_sl_data,
        training_prompt_data,
        training_sl_weights,
        training_prompt_weights,
    )
    efficiency_weights = reweighting.predicted_weights(
        bdt, test_prompt_data, test_prompt_weights
    )

    # Find prompt, SL, reweighted slices
    num_slices = 8
    bin_limits = (200, 1800)
    num_bins = 50
    plot_index = 1
    slice_index = 0

    Prompt_Slices = visualisations.Slices(
        num_slices, num_bins, bin_limits, plot_index, slice_index
    )
    SL_Slices = visualisations.Slices(
        num_slices, num_bins, bin_limits, plot_index, slice_index
    )
    Reweighted_Slices = visualisations.Slices(
        num_slices, num_bins, bin_limits, plot_index, slice_index
    )

    Prompt_Slices.add_points(test_prompt_data, test_prompt_weights)
    SL_Slices.add_points(test_sl_data, test_sl_weights)
    Reweighted_Slices.add_points(test_prompt_data, efficiency_weights)

    # Plot them
    visualisations.plot_slices(
        "python",
        (Prompt_Slices, SL_Slices, Reweighted_Slices),
        ("Prompt", "SL", "Reweighted"),
        ("red", "green", "blue"),
        r"M($\pi_1\pi_2$) /MeV"
    )


if __name__ == "__main__":
    #  plot_projections()
    #  n_calls = 250
    #  optimise(n_calls)
    # roc_score_test()
    plot_slices()
