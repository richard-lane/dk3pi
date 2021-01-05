import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import script_util

sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
import reweighting
import reweight_utils
import classification


def roc_score(prompt_points, sl_points, prompt_weights, sl_weights):
    """
    Compute ROC:AUC score

    """
    # Split data up
    points_train, points_test, labels_train, labels_test, weights_train, weights_test = classification.split_data_for_classification(
        prompt_points, sl_points, prompt_weights, sl_weights
    )

    # Train classifier
    classifier = classification.train_classifier(
        points_train, labels_train, weights_train
    )

    # Compute score
    return classification.classification_score(
        classifier, points_test, labels_test, weights_test
    )


def roc_score_test():
    """
    Compare ROC AUC scores for unweighted + reweighted distributions

    """
    # Read in data
    prompt_points, prompt_weights, sl_points, sl_weights = script_util.read_data()

    # Find the appropriate BDT weights
    training_prompt_data, test_prompt_data = np.array_split(prompt_points, 2)
    training_sl_data, test_sl_data = np.array_split(sl_points, 2)
    training_prompt_weights, test_prompt_weights = np.array_split(prompt_weights, 2)
    training_sl_weights, test_sl_weights = np.array_split(sl_weights, 2)
    print("Training BDT...")
    bdt = reweighting.init(
        training_sl_data,
        training_prompt_data,
        training_sl_weights,
        training_prompt_weights,
    )
    efficiency_weights = reweighting.predicted_weights(
        bdt, test_prompt_data, test_prompt_weights
    )

    linspace = np.linspace(0, 1)
    curve_before = classification.roc_curve(
        test_prompt_data, test_sl_data, test_prompt_weights, test_sl_weights
    )
    curve_after = classification.roc_curve(
        test_prompt_data, test_sl_data, efficiency_weights, test_sl_weights
    )
    plt.plot(curve_before[0], curve_before[1], label="Before Reweighting")
    plt.plot(curve_after[0], curve_after[1], label="After Reweighting")
    plt.plot(linspace, linspace, label="Indistinguishable", linestyle="--", color="k")
    plt.legend()
    plt.xlabel("FPR")
    plt.ylabel("TPR")
    plt.title("ROC Curves for Binary Classification of Prompt/Semileptonic Phsp Data")
    plt.savefig("roc_curve.png")
    plt.clf()


if __name__ == "__main__":
    roc_score_test()
