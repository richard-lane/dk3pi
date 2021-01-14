import matplotlib.pyplot as plt
import numpy as np
import sys
import os

sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
import reweight_utils


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


def read_data():
    """
    Returns phsp parametrised prompt_points, prompt_weights, sl_points, sl_weights

    """
    # Find phsp points for prompt + SL datasets
    print("Reading data...")
    prompt_points = reweight_utils.read_invariant_masses(
        "cut_wg_rs_prompt.root",
        "DecayTree",
        ("D0_P0_PX", "D0_P0_PY", "D0_P0_PZ", "D0_P0_PE"),
        ("D0_P1_PX", "D0_P1_PY", "D0_P1_PZ", "D0_P1_PE"),
        ("D0_P2_PX", "D0_P2_PY", "D0_P2_PZ", "D0_P2_PE"),
        ("D0_P3_PX", "D0_P3_PY", "D0_P3_PZ", "D0_P3_PE"),
    )
    sl_points = reweight_utils.read_invariant_masses(
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

    # Remove some SL points to make the distributions look more different
    # sl_indices_to_delete = []
    # for i in range(len(sl_points)):
    #    if sl_points[i][1] < 900 * np.random.random():
    #        sl_indices_to_delete.append(i)
    # sl_points = np.delete(sl_points, sl_indices_to_delete, axis=0)
    # sl_weights = np.delete(sl_weights, sl_indices_to_delete)

    # Remove some prompt points to test the BDT
    prompt_indices_to_delete = []
    for i in range(len(prompt_points)):
        if 800 < prompt_points[i][0] < 900:
            if 700 < prompt_points[i][1] < 800:
                if (
                    0.5 + (np.abs(prompt_points[i][1] - 750) / 125.0)
                    > np.random.random()
                ):
                    prompt_indices_to_delete.append(i)
    prompt_points = np.delete(prompt_points, prompt_indices_to_delete, axis=0)
    prompt_weights = np.delete(prompt_weights, prompt_indices_to_delete)

    return prompt_points, prompt_weights, sl_points, sl_weights


def hist_difference(hist1, hist2):
    """
    Takes a pair ofhistogram counts

    Returns (hist1 - hist2 counts), (hist1 - hist2 errors) assuming Poisson statistics

    """
    assert len(hist1) == len(hist2)

    # Find difference
    diff = np.subtract(hist1, hist2)

    # Find error on difference
    err1 = np.sqrt(hist1)
    err2 = np.sqrt(hist2)
    err = np.sqrt(np.add(err1 ** 2, err2 ** 2))

    return diff, err


def plot_hist_diffs(
    hist1, hist2, hist3, bin_centres, bin_widths, x_label, labels, diff_labels, title
):
    """
    Create a plot of three histograms, and a plot of the delta of hist1 and hist2 with hist3

    Errors assume Poisson statistics

    Doesn't show or save the plot- call plt.show() to show or plt.savefig() to save

    """
    fig, ax = plt.subplots(2, sharex=True, gridspec_kw={"height_ratios": [2.5, 1]})
    line_width = 0.5
    markersize = 0.0
    alpha = 0.5
    marker = "_"

    err1 = np.sqrt(hist1)
    err2 = np.sqrt(hist2)
    err3 = np.sqrt(hist3)

    # Plot the histograms
    for count, err, colour, label in zip(
        (hist1, hist2, hist3), (err1, err1, err3), ("red", "blue", "green"), labels
    ):
        ax[0].errorbar(
            bin_centres,
            count,
            yerr=err,
            xerr=bin_widths,
            label=label,
            alpha=alpha,
            color=colour,
            linewidth=line_width,
            marker=marker,
            markersize=markersize,
            fmt=" ",
        )

    # Find the differences between hists
    diff1, diff1_err = hist_difference(hist1, hist3)
    diff2, diff2_err = hist_difference(hist2, hist3)

    # Plot them
    for diff, err, colour, label in zip(
        (diff1, diff2), (diff1_err, diff2_err), ("red", "blue"), diff_labels
    ):
        ax[1].errorbar(
            bin_centres,
            diff,
            yerr=err,
            xerr=bin_widths,
            label=label,
            color=colour,
            linewidth=line_width,
            marker=marker,
            markersize=markersize,
            fmt=" ",
        )

    ax[0].legend()
    ax[1].legend()
    ax[0].set_ylabel("Counts")
    ax[1].set_ylabel(r"$\Delta$Counts")
    plt.xlabel(x_label)
    fig.subplots_adjust(hspace=0)
    fig.suptitle(title)
