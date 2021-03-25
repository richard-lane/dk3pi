"""
Utils for training, pickling and using bdts

"""
import hep_ml.reweight
import definitions
import pickle


def train(target, original, path, year, magnetisation, phsp_bin, decay_code):
    """
    Train a BDT reweighter on the provided data, and dump the pickled BDT to the provided path.

    :param target: target distribution (i.e. amplitude model)
    :param original: original distribution (i.e. LHCb MC)
    :param path: path to dump the pickled BDT to
    :param year: Data taking year, as an int
    :param magnetisation: magnetisation direction; "MagUp" or "MagDown"
    :param phsp_bin: phsp bin number according to definitions.PHSP_BINS. Indexed from 0
    :param decay_code: 27165070, 27165071 or 27165072

    """
    # Train the BDT
    bdt = hep_ml.reweight.GBReweighter(n_estimators=200)
    bdt.fit(original=original, target=target)

    # Pickle it
    with open(path, "wb") as f:
        pickle.dump(bdt, f)

    # Save the location of the pickled BDT to the dict of BDTs
    with open(definitions.BDT_PATHS, "rb") as f:
        paths = pickle.load(f)

    paths[year][magnetisation][phsp_bin] = path
    with open(definitions.BDT_PATHS, "wb") as f:
        pickle.dump(paths, f)


def find_weights(bdt_path, data, weights=None):
    """
    Find the weights to apply to data using the bdt pickled at bdt_path

    NB: the weights are not normalised

    :param bdt_path: location of file holding a trained BDT. Will be unserialised using pickle.load
    :param data: the data to reweight
    :param weights: weights of data before reweighting (e.g. sWeights)

    """
    with open(bdt_path, "rb") as f:
        bdt = pickle.load(f)

    return bdt.predict_weights(data, original_weight=weights)
