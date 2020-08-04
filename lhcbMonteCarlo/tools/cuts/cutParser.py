"""
Script for reading a module of cuts implemented in python and applying them to a dataset

Needs to import a module from the command line, find what cuts it contains and apply them to our dataset

Assumes cut functions start cut*

"""
import argparse
import importlib
import inspect
import uproot


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
    cut_names = tuple(cut for cut in dir(cuts_lib) if cut.startswith("cut"))

    # For these cuts, find the string representation of their arguments
    cut_args = tuple(
        tuple(inspect.getfullargspec(getattr(cuts_lib, cut)).args) for cut in cut_names
    )

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
    print(len(data["Dstar_M"]))


if __name__ == "__main__":
    main(cli())
