"""
Stuff that's helpful for reweighting using python BDT

Python I/O, phase space parametrisation and BDT hyperparameter optimisation

"""
import uproot


def read_branch(file_name: str, tree_name: str, branch_name: str) -> np.ndarray:
    """
    Read the contents of a ROOT branch into a numpy.ndarray

    """
    tree = uproot.open(file_name)[tree_name]

    return tree.array(branch_name)
