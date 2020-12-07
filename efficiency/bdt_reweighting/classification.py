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


def roc_curve(prompt_points, sl_points, prompt_weights, sl_weights):
    """
    Return FPR, TPR, threshholds suitable for making a ROC curve plot

    """
    points_train, points_test, labels_train, labels_test, weights_train, weights_test = split_data_for_classification(
        prompt_points, sl_points, prompt_weights, sl_weights
    )

    # Train classifier
    classifier = train_classifier(points_train, labels_train, weights_train)

    # Find probabilities for the test data being correctly classified
    probs = classifier.predict_proba(points_test)[:, 1]

    # Find which weights correspond to prompt/sl
    s_weights = weights_test * (labels_test == 1)
    p_weights = weights_test * (labels_test == 0)

    # Find probabilities of false classification and decision function values (i think)
    threshhold, probs = np.unique(probs, return_inverse=True)

    # Find cumulative type 1/2 error rates
    tpr = np.bincount(probs, weights=s_weights)[::-1].cumsum()
    fpr = np.bincount(probs, weights=p_weights)[::-1].cumsum()

    # Normalise probabilities
    tpr /= tpr[-1]
    fpr /= fpr[-1]

    return fpr, tpr, threshhold[::-1]
