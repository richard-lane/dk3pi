"""
Stuff for quantifying and visualising how well our BDT reweighting is going

"""
import numpy as np
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split


def _split_data_for_classification(
    original_distribution, target_distribution, original_weights, target_weights
):
    """
    Combine the original + target distributions into a single sample with labels, and then split using sklearn.model_selection.train_test_split

    Uses 0 to label original distribution and 1 for the target distribution

    """
    assert len(original_distribution) == len(original_weights)
    assert len(target_distribution) == len(target_weights)

    data = np.concatenate([original_distribution, target_distribution])
    weights = np.concatenate(original_weights, target_weights)

    # Use 0 to label the original distribution and 1 to label the target
    data = np.concatenate([original_distribution, target_distribution])
    labels = np.array([0] * len(original_distribution) + [1] * len(target_distribution))

    return train_test_split(data, labels, weights)


def train_classifier(distribution, classification, weights):
    """
    Train a classifier on a distribution of points + a classification of those points

    original_weights should attempt to take original_distribution to target_distribution * target_weights

    Returns whatever sklearn.ensemble.GradientBoostingClassifier.fit returns

    """
    return GradientBoostingClassifier().fit(distribution, classification, weights)
