"""
Script for reading a module of cuts implemented in python and applying them to a dataset

Needs to import a module from the command line, find what cuts it contains and apply them to our dataset

Assumes cut functions start cut*

"""
import argparse
import importlib
import inspect
import numpy as np
import uproot


def write_root_file(data, filename, tree_name):
    """
    Write a python dict of data to a ROOT file

    """
    # Create a dict of {branchname : type} that we use to declare our tree
    data_types = {key: type(data[key][0]) for key in data}
    tree = uproot.newtree(data_types)

    # Declare the ROOT file
    print(f"Creating ROOT file '{filename}'")
    root_file = uproot.recreate(filename)
    root_file[tree_name] = tree

    # Fill tree
    print(f"Populating ROOT file with tree '{tree_name}'")
    root_file[tree_name].extend(data)


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
    tree = uproot.open(input_file)[decay_tree]

    # If branches not provided then just use all of them
    read_branches = branches
    if not branches:
        read_branches = [branch.decode("utf-8") for branch in tree.allkeys()]

    data = tree.arrays(read_branches, namedecode="utf-8")

    # Uproot currently only supports writing using
    #     int8, float64, float32, int32, int64, bool or int16
    # i.e. cannot handle unsigned things
    for key in data:
        if np.dtype(type(data[key][0])).kind == "u":
            raise NotImplementedError(
                f"\nWill be unable to write branch {key} of type {type(data[key][0])}\n"
                "Writing unsigned integers to file currently unsupported by uproot.\n"
                "I might add a feature to cast unsigned->signed branches if needed"
            )

    return data


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
        "--out_file",
        help="File to write to, defaults to 'out.root'",
        default="out.root",
    )
    parser.add_argument(
        "--decay_tree",
        help="Decay Tree object to read, defaults to TupleDstToD0pi_D0ToKpipipi/DecayTree for no reason",
        default="TupleDstToD0pi_D0ToKpipipi/DecayTree",
    )
    parser.add_argument(
        "--out_tree",
        help="Tree to write to; Uproot does not currently support subdirectories so this can't be nested",
        default="myTree",
    )

    return parser.parse_args()


def main(args):
    # Find our cut functions and the arguments they take
    cuts_lib = importlib.import_module(args.cuts_module)
    cut_fcns, cut_args = cut_functions(cuts_lib)

    # If no branch name is provided in the config module, set the branches variable to None
    # The read function will know to just read all branches
    try:
        branches = cuts_lib.BRANCHES
    except AttributeError:
        branches = None
    data = read_root_branches(args.in_file, args.decay_tree, branches)

    # Perform the cuts and write to our output file
    perform_cuts(data, cut_args, cut_fcns)
    write_root_file(data, args.out_file, args.out_tree)


if __name__ == "__main__":
    main(cli())
