"""
Illustrative graphs of BDT reweighting

"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# This is really horrible but i can't think of a better way of doing it
# Ideally i'd like to set a global python include path via the CMake build system...
sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
import reweighting


def find_weights(generated_events, training_detected_events, detected_events):
    """
    Find the weights to apply to detected_events given a target distribution generated_events and a measured distribution training_detected_events

    """
    # Train a BDT on the detected + measured events
    print("Training BDT...")
    bdt = reweighting.init(generated_events, training_detected_events)

    # We want to reweight detected_events to look like generated_events
    # => We want to have len(generated_events) events in our reweigted distribution
    num_events = len(generated_events)

    # Given another measured sample, find the weights that apply to it
    print("finding weights")
    return reweighting.predicted_weights(bdt, detected_events, num_events)


def reweight_1d_gaussian(num_events):
    """
    Create plots for reweighting one 1d gaussian to another

    """
    # Create a target dataset
    target_mean = 5.0
    measured_mean = 6.0
    std_dev = 1.0
    print("Generating target data")
    generated_events = np.random.normal(target_mean, std_dev, size=(num_events, 1))

    # Create two observed datasets- one for training, one to reweight
    print("Generating measured data")
    num_detected_events = num_events // 2
    detected_events = np.random.normal(
        measured_mean, std_dev, size=(num_detected_events, 1)
    )
    more_detected_events = np.random.normal(
        measured_mean, std_dev, size=(num_detected_events, 1)
    )

    weights = find_weights(generated_events, detected_events, more_detected_events)

    # Plot datasets, including the reweighted one
    num_bins = 100
    bins = np.linspace(0, 10, num_bins + 1)
    plt.hist(generated_events, bins=bins, histtype="step")
    plt.hist(more_detected_events, bins=bins, histtype="step")
    plt.hist(more_detected_events, bins=bins, weights=weights, histtype="step")
    plt.legend(("MC", "Detected", "Reconstructed"))
    plt.show()


def reweight_multidimensional_gaussian(num_events):
    """
    Create plots for reweighting a multidimensional gaussian

    """
    # diagonal covariance matrix => uncorrelated
    covariance_matrix = np.array([[2.0, -1.0, 0.5], [-1.0, 1.0, -0.6], [0.5, -0.6, 1.5]])
    target_means = np.array([3, 5.4, 4.6])
    measured_means = np.array([2, 6, 4.6])
    dimensionality = len(covariance_matrix)

    # Generate multidimensional Gaussian for target data
    print("Generating target data")
    generated_events = np.random.multivariate_normal(
        target_means, covariance_matrix, size=num_events, check_valid="raise"
    )

    # Generate two multidimensional Gaussians for training and observed data
    print("Generating measured data")
    num_detected_events = num_events // 2
    detected_events = np.random.multivariate_normal(
        measured_means, covariance_matrix, size=num_detected_events, check_valid="raise"
    )
    more_detected_events = np.random.multivariate_normal(
        measured_means, covariance_matrix, size=num_detected_events, check_valid="raise"
    )

    # Train a BDT on the detected + measured events
    weights = find_weights(generated_events, detected_events, more_detected_events)

    # Plot projections of target, training_measured + reweighted_measured data
    num_bins = 100
    bins = np.linspace(-5, 10, num_bins + 1)
    fig = plt.figure()
    fig.suptitle("3d Gaussian Projections")
    for i in range(3):
        plt.subplot(1, 3, i + 1)
        plt.hist(generated_events[:, i], bins=bins, histtype="step")
        plt.hist(more_detected_events[:, i], bins=bins, histtype="step")
        plt.hist(
            more_detected_events[:, i], bins=bins, weights=weights, histtype="step"
        )
        plt.title(f"Projection {i}")

    fig.legend(("MC", "Detected", "Reconstructed"), loc="center right")
    plt.show()


def main():
    # reweight_1d_gaussian(10 ** 5)
    reweight_multidimensional_gaussian(10 ** 5)


if __name__ == "__main__":
    main()
