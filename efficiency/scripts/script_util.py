import numpy as np
import sys
import os

sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
import reweight_utils


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

    # Remove some points to make the distributions look more different
    indices_to_delete = []
    for i in range(len(sl_points)):
        if sl_points[i][1] < 900 * np.random.random():
            indices_to_delete.append(i)

        # Remove some more points to test the BDT
        elif 600 < sl_points[i][1] < 800:
            #    if 800 < sl_points[i][0] < 900:
            if (0.2 + np.abs(sl_points[i][1] - 700) / 125.0) > np.random.random():
                indices_to_delete.append(i)

    sl_points = np.delete(sl_points, indices_to_delete, axis=0)
    sl_weights = np.delete(sl_weights, indices_to_delete)

    return prompt_points, prompt_weights, sl_points, sl_weights
