"""
Stuff for quantifying and visualising how well our BDT reweighting is going

"""
import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split


def split_data_for_classification(
    original_distribution, target_distribution, original_weights, target_weights
):
    """
    Combine the original + target distributions into a single sample with labels, and then split using sklearn.model_selection.train_test_split

    Uses 0 to label original distribution and 1 for the target distribution

    """
    assert len(original_distribution) == len(original_weights)
    assert len(target_distribution) == len(target_weights)

    data = np.concatenate([original_distribution, target_distribution])
    weights = np.concatenate((original_weights, target_weights))

    # Use 0 to label the original distribution and 1 to label the target
    labels = np.array([0] * len(original_distribution) + [1] * len(target_distribution))

    return train_test_split(data, labels, weights)


def train_classifier(distribution, classification, weights):
    """
    Train a classifier on a distribution of points + a classification of those points

    original_weights should attempt to take original_distribution to target_distribution * target_weights

    Returns whatever sklearn.ensemble.GradientBoostingClassifier.fit returns

    """
    return GradientBoostingClassifier().fit(distribution, classification, weights)


def classification_score(classifier, data, classifications, weights):
    """
    Find the ROC AUC (Area Under Curve: Receiver Operating Characteristic) score given a trained classifier, a dataset, its classification and weights

    Data probably has to be in some sort of order

    """
    return roc_auc_score(
        classifications, classifier.predict_proba(data)[:, 1], sample_weight=weights
    )
