"""
Tune the efficiency BDT using a combined chi squared score

"""
import numpy as np
import matplotlib.pyplot as plt
import skopt
import sys
import os
import warnings
import script_util

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


def optimise(n_calls):
    """
    Find the optimal BDT parameters for a given dataset...

    Try n_calls BDT configurations; starts with n_calls//5 random points

    """
    # Read data from files + perform phsp parametrisation
    prompt_points, prompt_weights, sl_points, sl_weights = script_util.read_data()

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

if __name__ == "__main__":
    n_calls = 10
    optimise(n_calls)