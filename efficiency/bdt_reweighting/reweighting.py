import numpy as np
from hep_ml.reweight import GBReweighter


def init(
    mc_data,
    real_data,
    mc_weights=None,
    real_data_weights=None,
    n_estimators=40,
    learning_rate=0.2,
    max_depth=3,
    min_samples_leaf=200,
    loss_regularization=5.0,
    gb_args=None,
):
    """
    Initialise, train + return a BDT for performing the D->K3pi efficiency estimate

    Pass in an array-like of data points, each of which is an array-like of co-ordinates
    e.g. [[1, 2, 3, 4], [1, 4, 2, 5], ...]

    Also optionally pass in an array-like of weights; if unweighted pass in None

    """
    reweighter = GBReweighter(
        n_estimators,
        learning_rate,
        max_depth,
        min_samples_leaf,
        loss_regularization,
        gb_args,
    )
    reweighter.fit(
        original=real_data,
        target=mc_data,
        original_weight=real_data_weights,
        target_weight=mc_weights,
    )

    return reweighter


def predicted_weights(
    reweighter: GBReweighter, points, weights=None, expected_num_points=None
):
    """
    Find the predicted weights at a given point from a (trained) reweighter

    If we want our reweighted distribution to have a certain number of points, provide this above
    This bit doesn't actually work right so just scale it yourself later

    """
    eff_weights = reweighter.predict_weights(points, weights)

    # We may want to fiddle with the weights so we get the right normalisation afterwards
    normalisation = 1.0
    if expected_num_points:
        avg_weight = np.mean(eff_weights)
        normalisation *= expected_num_points / (len(points) * avg_weight)

    return normalisation * eff_weights
