"""
BDT reweighting in pure python

"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import warnings

# This is really horrible but i can't think of a better way of doing it
# Ideally i'd like to set a global python include path via the CMake build system...
sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
import reweighting
import reweight_utils


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


def main():
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

    # Compare the reweighted test prompt + test SL data by plotting some histograms
    bins = np.linspace(200, 1800, 250)
    for i in range(5):
        prompt, edges = np.histogram(
            test_prompt_data[:, i], bins=bins, weights=test_prompt_weights
        )
        reweighted, _ = np.histogram(
            test_prompt_data[:, i], bins=bins, weights=efficiency_weights
        )
        sl, _ = np.histogram(test_sl_data[:, i], bins=bins, weights=test_sl_weights)

        # Find errors
        prompt_err = np.sqrt(prompt)
        reweighted_err = np.sqrt(reweighted)
        sl_err = np.sqrt(sl)

        # Rescale histograms and errors
        prompt_integral = np.sum(prompt)
        prompt /= prompt_integral
        prompt_err /= prompt_integral

        sl_integral = np.sum(sl)
        sl /= sl_integral
        sl_err /= sl_integral

        reweighted_integral = np.sum(reweighted)
        reweighted /= reweighted_integral
        reweighted_err /= reweighted_integral

        # Find bin centres
        centres = np.mean(np.vstack([edges[0:-1], edges[1:]]), axis=0)

        # Make plots
        plt.errorbar(
            centres,
            prompt,
            yerr=prompt_err,
            label="Prompt",
            color="red",
            linewidth=0.5,
            marker=".",
            markersize=0.5,
        )
        plt.errorbar(
            centres,
            reweighted,
            yerr=reweighted_err,
            label="Reweighted",
            color="blue",
            linewidth=0.5,
            marker=".",
            markersize=0.5,
        )
        plt.errorbar(
            centres,
            sl,
            yerr=sl_err,
            label="SL",
            color="green",
            marker=".",
            linewidth=0.5,
            markersize=0.5,
        )
        plt.legend()

        plt.savefig(f"{i}.png", format="png", dpi=1000)
        plt.clf()

        # Print chi squareds
        print(f"Orig: {chi_squared(prompt, sl)}\tReweighted: {chi_squared(reweighted, sl)}")


if __name__ == "__main__":
    main()
