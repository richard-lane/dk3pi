"""
BDT reweighting in pure python

"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import script_util


# This is really horrible but i can't think of a better way of doing it
# Ideally i'd like to set a global python include path via the CMake build system...
sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
import reweighting
import reweight_utils
import visualisations


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
    # for count, err, colour, label in zip(
    #    (prompt_counts, reweighted_counts, sl_counts),
    #    (prompt_err, reweighted_err, sl_err),
    #    ("red", "blue", "green"),
    #    ("Prompt", "Reweighted", "SL"),
    # ):
    for count, err, colour, label in zip(
        (prompt_counts, sl_counts),
        (prompt_err, sl_err),
        ("red", "green"),
        ("Prompt", "SL"),
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
        prompt, sl, reweighted = script_util.bin_data(
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
        script_util.rescale(prompt, prompt_err)
        script_util.rescale(sl, sl_err)
        script_util.rescale(reweighted, reweighted_err)

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


def read_and_reweight():
    """
    Read data from ROOT files, train BDT + find weights

    Trains the BDT on half the data; returns the other half

    Returns prompt points, prompt weights, sl points, sl_weights, efficiency weights

    """
    # Read data from files + perform phsp parametrisation
    prompt_points, prompt_weights, sl_points, sl_weights = script_util.read_data()

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

    return (
        test_prompt_data,
        test_prompt_weights,
        test_sl_data,
        test_sl_weights,
        efficiency_weights,
    )


def plot_projections(
    prompt_points, prompt_weights, sl_points, sl_weights, efficiency_weights
):
    """
    Plot phsp projections with default BDT

    """

    # Plot phsp projections
    bins = np.linspace(200, 1800, 250)
    make_plots(
        prompt_points, sl_points, prompt_weights, sl_weights, efficiency_weights, bins
    )


def plot_slices(
    prompt_points, prompt_weights, sl_points, sl_weights, efficiency_weights
):
    # Find prompt, SL, reweighted slices
    num_slices = 12
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

    Prompt_Slices.add_points(prompt_points, prompt_weights)
    SL_Slices.add_points(sl_points, sl_weights)
    Reweighted_Slices.add_points(prompt_points, efficiency_weights)

    # Plot them
    visualisations.plot_slices(
        "python",
        (Prompt_Slices, SL_Slices, Reweighted_Slices),
        ("Prompt", "SL", "Reweighted"),
        ("red", "green", "blue"),
        r"M($\pi_1\pi_2$) /MeV",
        r"M($K\pi_1$) /MeV",
    )


if __name__ == "__main__":
    prompt_points, prompt_weights, sl_points, sl_weights, efficiency_weights = (
        read_and_reweight()
    )
    plot_projections(
        prompt_points, prompt_weights, sl_points, sl_weights, efficiency_weights
    )
    plot_slices(
        prompt_points, prompt_weights, sl_points, sl_weights, efficiency_weights
    )

