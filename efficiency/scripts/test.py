"""
Train a BDT, pickle it, unpickle it, use it

"""
import hep_ml.reweight
import numpy as np
import matplotlib.pyplot as plt
import pickle


def _generate():
    N = 100000
    print("gen target")
    target = np.random.multivariate_normal((0.0, 0.0), ((1.0, 0.0), (0.0, 1.0)), N)

    # Generate a slightly different 2d gaussian
    print("gen original")
    original = np.random.multivariate_normal((0.3, 0.0), ((1.1, 0.1), (0.1, 0.9)), N)
    return target, original


def _train_bdt():
    """
    Returns a trained BDT

    """
    target, original = _generate()

    # Train a BDT to reweight
    print("train bdt")
    bdt = hep_ml.reweight.GBReweighter()
    bdt.fit(original=original, target=target)

    return bdt


def _apply_pickled_bdt():
    # Unpickle bdt
    with open("test_bdt.pickle", "rb") as f:
        bdt = pickle.load(f)

    # Use it to reweight this data
    target, original = _generate()
    weights = bdt.predict_weights(original)

    kw = {"bins": 100, "alpha": 0.3, "density": True}
    for i in range(2):
        plt.hist(target[:, i], **kw)
        plt.hist(original[:, i], **kw)
        plt.show()

    for i in range(2):
        plt.hist(target[:, i], **kw)
        plt.hist(original[:, i], weights=weights, **kw)
        plt.show()


def main():
    # Train a bdt
    bdt = _train_bdt()

    # Pickle it
    print("pickling bdt")
    with open("test_bdt.pickle", "wb") as f:
        pickle.dump(bdt, f)

    # Do something else
    print("the bdt has been pickled")

    # Use it to reweight
    _apply_pickled_bdt()


if __name__ == "__main__":
    main()
