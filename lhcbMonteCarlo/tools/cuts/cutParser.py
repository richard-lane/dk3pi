"""
Script for reading a module of cuts implemented in python and applying them to a dataset

Needs to import a module from the command line, find what cuts it contains and apply them to our dataset

Assumes cut functions start cut*

"""
import argparse
import importlib
import inspect
import uproot
import numpy as np


def perform_cuts(data, cut_branches, cuts):
    """
    Perform the specified cuts on a dataset

    Data should be a dict of {branchname: data}

    branches should be an iterable of tuples that the cuts apply to

    e.g. to perform a cut that needs to know about Dstar and D0 mass, call something like
        perform_cuts(data, [("Dstar_M", "D0_M")], [my_cut_fcn])

    """
    print(f"Performing cuts using:")
    for branch_tuple, function in zip(cut_branches, cuts):
        print("\t" + f"{function.__name__}{branch_tuple}")

    # Delete data from our datasets according to our cuts
    num_datapoints = len(data[cut_branches[0][0]])
    num_cuts = len(cuts)

    # This isn't optimised
    # Need to have an ordered list of our indices, as we will want to remove in the order highest->lowest when it comes to actually performing the cuts
    indices_to_remove = []
    for i in range(num_datapoints):
        for cut in range(num_cuts):
            cut_args = tuple(data[arg][i] for arg in cut_branches[cut])
            if i not in indices_to_remove:
                if not cuts[cut](*cut_args):
                    indices_to_remove.append(i)

    print(f"Removing {len(indices_to_remove)} events of {num_datapoints}...")

    for i in data:
        # Iterate through our indices in descending order so we don't accidentally skip any values
        for index in np.flip(indices_to_remove):
            # This is likely HUGELY slow, since numpy arrays are immutable and so we copy the entire array after each delete probably
            data[i] = np.delete(data[i], index)

    print(f"After cuts: {len(data[cut_branches[0][0]])} data points")


def read_root_branches(input_file, decay_tree, branches):
    """
    Read the provided ROOT branches into numpy arrays or something

    Returns a dict of arrays i guess

    """
    root_file = uproot.open(input_file)
    return root_file[decay_tree].arrays(branches, namedecode="utf-8")


def cut_functions(cuts_lib):
    """
    From our cuts module, identify the functions and arguments needed to perform cuts

    Returns a tuple of cut functions and tuple of tuples of their arguments

    """

    # Find the functions in our module that start with 'cut'
    cut_names = tuple(
        getattr(cuts_lib, cut) for cut in dir(cuts_lib) if cut.startswith("cut")
    )

    # For these cuts, find the string representation of their arguments
    cut_args = tuple(tuple(inspect.getfullargspec(cut).args) for cut in cut_names)

    return cut_names, cut_args


def cli():
    """
    Argument parsing; returns an argparse ArgumentParser object

    """
    parser = argparse.ArgumentParser(
        description="Cutting unwanted data from a ROOT file. Cuts provided by a user-provided cuts module"
    )
    parser.add_argument("in_file", help="ROOT file to read")
    parser.add_argument("cuts_module", help="Python module defining cuts")
    parser.add_argument(
        "--outFile", help="File to write to, defaults to 'out.root'", default="out.root"
    )
    parser.add_argument(
        "--decay_tree",
        help="Decay Tree object to read, defaults to TupleDstToD0pi_D0ToKpipipi/DecayTree for no reason",
        default="TupleDstToD0pi_D0ToKpipipi/DecayTree",
    )

    return parser.parse_args()


def main(args):
    cuts_lib = importlib.import_module(args.cuts_module)
    cut_fcns, cut_args = cut_functions(cuts_lib)

    branches = cuts_lib.BRANCHES

    data = read_root_branches(args.in_file, args.decay_tree, branches)

    perform_cuts(data, cut_args, cut_fcns)


if __name__ == "__main__":
    main(cli())
