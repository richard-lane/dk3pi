import matplotlib.pyplot as plt
from scipy.stats import chi2 as scipy_chi2
import numpy as np
import sys
import os
import phasespace
import warnings
import tqdm

sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
import reweight_utils


def flat_phsp_points(num_decays):
    """
    Return the K, pi1, pi2, pi3 kinematic parameters for a number of D->K3pi decays, uniformly distributed in phase space

    Does not order the pions in any way; pi1/2/3 are indistinguishable

    """
    pi_mass = 139.570
    k_mass = 493.677
    d_mass = 1864.84
    generator = phasespace.nbody_decay(
        d_mass, (k_mass, pi_mass, pi_mass, pi_mass), names=("K", "pi1", "pi2", "pi3")
    )

    # Initialise particle arrays to the right shape
    k = np.zeros((4, num_decays))
    pi1 = np.zeros((4, num_decays))
    pi2 = np.zeros((4, num_decays))
    pi3 = np.zeros((4, num_decays))

    # Seems like the maximum weight it is possible to generate is <0.12
    max_weight = 0.12

    # Progress bar
    with tqdm.tqdm(total=num_decays) as pbar:
        num_generated = 0

        while num_generated < num_decays:
            # Generate a chunk of our desired number of decays, since generating many particles is almost as fast as generating 1
            chunk_size = num_decays * 2
            weights, particles = generator.generate(chunk_size, normalize_weights=True)

            # Generate a load of random numbers to compare to
            random_numbers = max_weight * np.random.random(chunk_size)

            # Iterate over the particles and weights we generated, exiting early if we overflow
            for rnd, weight, this_k, this_pi1, this_pi2, this_pi3 in zip(
                random_numbers,
                weights,
                particles["K"].numpy(),
                particles["pi1"].numpy(),
                particles["pi2"].numpy(),
                particles["pi3"].numpy(),
            ):
                # Don't want to cut off our distribution by accident
                assert weight < max_weight

                # Check if our point is accepted and if so insert it into the array
                if rnd < weight:
                    try:
                        k.T[num_generated] = this_k
                        pi1.T[num_generated] = this_pi1
                        pi2.T[num_generated] = this_pi2
                        pi3.T[num_generated] = this_pi3

                        num_generated += 1
                        pbar.update(1)

                    # If we overshoot the end of array, stop
                    except IndexError:
                        break

                    # Would maybe be faster to generate smaller and smaller chunks as we iterate, but this is fine

    return k, pi1, pi2, pi3


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


def _bin_contents_and_errs(data, weights, bins):
    """
    Helper for chisq fcn

    Returns weighted bin contents and the contribution to the error term

    i.e. returns (bin contents, bin errors) where bin errors are the sum of the (weights squared)

    Raises ValueError in case of under or overflows

    """
    num_bins = len(bins) - 1
    bin_contents = np.zeros(num_bins)
    bin_errs = np.zeros(num_bins)

    # This won't be as fast as using numpy, but will be fine for my purposes
    for x_i, w_x in zip(data, weights):
        # Find which bin this point belongs in
        bin_index = np.digitize(x_i, bins) - 1
        if bin_index == -1 or bin_index == num_bins:
            raise ValueError(
                f"Bin over/underflow with point {x_i} and bin extrema {bins[0]},{bins[-1]}"
            )

        # Increment the right bin content
        bin_contents[bin_index] += w_x

        # Add weight squared to the error term for this bin
        bin_errs[bin_index] += w_x ** 2

    return bin_contents, bin_errs


def chi_sq_distance(x, y, bins, x_weights=None, y_weights=None):
    """
    Find the chi squared distance between two 1d distributions x and y.

    Bin the distributions x and y into the provided bins
    and calculate the chi squared distance between them using the following formula:
        $\chi^2 = \sum_i \frac{(x_i - y_i)^2}{\sum_{bin}w_x^2 + \sum_{bin}w_y^2}$
    Where e.g. $w_x$ are x's weights in the given bin.
    Warns if there are fewer than 10 (weighted) events in any bin.

    :param x: a 1d iterable of data values in a distribution.
    :param y: a 1d iterable of data values in a distribution.
    :param bins: the bins to insert x and y into. For N bins, should have N+1 entries: i.e. must contain every bin edge.
    :param x_weights: the weights to apply to x.
    :param y_weights: the weights to apply to y.

    :returns: the value of chisq calculated. Not normalised
    :returns: the corresponding p value

    :raises: `ValueError` if len(x) != len(x_weights) or len(y) != len(y_weights)
    :raises: `ValueError` if a points under or overflows the provided bins

    """
    if x_weights is not None and len(x) != len(x_weights):
        raise ValueError(
            f"x data and weights have different lengths: {len(x)} and {len(x_weights)}"
        )

    if y_weights is not None and len(y) != len(y_weights):
        raise ValueError(
            f"y data and weights have different lengths: {len(y)} and {len(y_weights)}"
        )

    # Possibly slightly (memory) wasteful
    if x_weights is None:
        x_weights = np.ones_like(x)
    if y_weights is None:
        y_weights = np.ones_like(y)

    # Find the contents and errors associated with each bin
    x_contents, x_errs = _bin_contents_and_errs(x, x_weights, bins)
    y_contents, y_errs = _bin_contents_and_errs(y, y_weights, bins)

    # If any bin contains fewer than 10 (weighted) points, throw a warning
    num_bins = len(bins) - 1
    for contents, label in zip((x_contents, y_contents), ("x", "y")):
        for content, i in zip(contents, range(num_bins)):
            if content < 10:
                warnings.warn(
                    f"Fewer than 10 points ({content}) found in bin {i} for {label} distribution; chisq distance may be unreliable",
                    category=RuntimeWarning,
                )

    # Sum up the bin differences and divide by bin errors to find chi squared
    bin_diffs_sq = np.square(np.subtract(y_contents, x_contents))
    bin_errs = np.add(x_errs, y_errs)
    chisq = np.sum(np.divide(bin_diffs_sq, bin_errs))

    # We want the two-tailed p-value ? TODO
    p_value = scipy_chi2.sf(chisq, num_bins)

    return chisq, p_value

