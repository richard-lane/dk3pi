import numpy as np
from hep_ml.reweight import GBReweighter


def init(mc_data, real_data, mc_weights=None, real_data_weights=None):
    """
    Initialise, train + return a BDT for performing the D->K3pi efficiency estimate

    Pass in an array-like of data points, each of which is an array-like of co-ordinates
    e.g. [[1, 2, 3, 4], [1, 4, 2, 5], ...]

    Also optionally pass in an array-like of weights; if unweighted pass in None

    """
    # reweighter = GBReweighter(max_depth=7, gb_args={"max_features": 1})
    reweighter = GBReweighter()
    reweighter.fit(
        original=real_data,
        target=mc_data,
        original_weight=real_data_weights,
        target_weight=mc_weights,
    )

    return reweighter


def predicted_weights(reweighter: GBReweighter, points, expected_num_points=None):
    """
    Find the predicted weights at a given point from a (trained) reweighter

    If we want our reweighted distribution to have a certain number of points, provide this above

    """
    weights = reweighter.predict_weights(points)

    # We may want to fiddle with the weights so we get the right normalisation afterwards
    normalisation = 1.0
    if expected_num_points:
        avg_weight = np.mean(weights)
        normalisation *= expected_num_points / (len(points) * avg_weight)

    return normalisation * weights
