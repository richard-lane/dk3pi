"""
Library of functions for calculating and plotting a thing

"""
from numpy import concatenate, linspace, histogram
import sys, os

sys.path.append(
    os.path.abspath(
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "bdt_reweighting")
    )
)
import script_util


def goodness_of_fit(
    original_distribution,
    target_distribution,
    Classifier,
    bins,
    original_weights=None,
    target_weights=None,
):
    """
    Analyze the multidimensional goodness of fit between two distributions

    Takes a (trained) binary classifier and classifies the original + target distributions.
    Returns the probability of classifying the original/target distributions as TODO as what??

    :param original_distribution: a numpy.ndarrary of shape (N, d) for N d-dimensional points in the original distribution
    :param target_distribution: a numpy.ndarray of shape (M, d) for M d-dimensional points in the target distribution (i.e. the distribution to compare original_distribution to)
    :param Classifier: a binary classifer that has been trained to distinguish the two distributions. Must provide a predict_proba method
    :param bins: bin limits to use when calculating chi squared distance. Should contain the left bin edge of each bin, plus the rightmost edge of the highest bin.
    :param original_weights: weights to apply to the original distribution's probabilities when calculating chi squared
    :param target_weights: weights to apply to the target distribution's probabilities when calculating chi squared

    :returns: classification probability of the original distribution
    :returns: classification probability of the target distribution
    :returns: chi squared between these probabilities (based on hard-coded bins)
    :returns: p-value

    """
    # Find the classification probabiltiies
    orig_prob, target_prob = (
        Classifier.predict_proba(original_distribution)[:, 0],
        Classifier.predict_proba(target_distribution)[:, 0],
    )

    # find a chi squared and p value between these two histograms
    chisq, p = script_util.chi_sq_distance(
        orig_prob, target_prob, bins, original_weights, target_weights
    )

    return orig_prob, target_prob, chisq, p
