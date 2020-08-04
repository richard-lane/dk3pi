"""
Script for parsing our LHCb Monte Carlo dataset and just extracting the data we want/need

Implments the cuts in cuts.txt; uses uproot

"""
import argparse
import numpy as np
import uproot4
import uproot


def d_mass_cut(d_masses, i):
    """
    Returns bool for our D mass cut
    True if d_masses[i] is in the allowed range

    d_masses in MeV

    """
    return 1840.83 < d_masses[i] < 1888.83


def delta_m_cut(d_star_masses, d_masses, i):
    """
    Returns bool for our DeltaM cut
    True if deltam[i] is in the allowed range

    masses in MeV

    """
    return 139.3 < d_star_masses[i] - d_masses[i] < 152


def d_impact_param_cut(ip_chisqs, i):
    """
    Returns bool for impact parameter cut
    True if our chisq is below the allowed value

    Takes in D0_IPCHI2_OWNPV

    """
    return ip_chisqs[i] < 9


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


def perform_cuts(data, branches, cuts):
    """
    Perform the specified cuts on a dataset

    Data should be a dict of {branchname: data}

    branches should be an iterable of tuples that the cuts apply to

    e.g. to perform a cut that needs to know about Dstar and D0 mass, call something like
        perform_cuts(data, [("Dstar_M", "D0_M")], [my_cut_fcn])

    """
    if len(branches) != len(cuts):
        raise Exception("Num branches != num cuts provided")

    print(f"Performing cuts using variables:")
    for branch_tuple in branches:
        print("\t" + f"{branch_tuple}")

    # Delete data from our datasets according to our cuts
    num_datapoints = len(data[branches[0][0]])
    num_cuts = len(cuts)

    # Unnecessary copy but shouldn't make too much difference
    branch_data = [
        [data[branch] for branch in branch_tuple] for branch_tuple in branches
    ]

    # This isn't optimised
    # Need to have an ordered list of our indices, as we will want to remove in the order highest->lowest when it comes to actually performing the cuts
    indices_to_remove = []
    for i in range(num_datapoints):
        for cut in range(num_cuts):
            if i not in indices_to_remove:
                # Could probably use np.ravel
                if not cuts[cut](*(branch_data[cut]), i):
                    indices_to_remove.append(i)

    print(f"Removing {len(indices_to_remove)} events of {num_datapoints}...")

    for i in data:
        # Iterate through our indices in descending order so we don't accidentally skip any values
        for index in np.flip(indices_to_remove):
            # This is likely HUGELY slow, since numpy arrays are immutable and so we copy the entire array after each delete
            data[i] = np.delete(data[i], index)

    print(f"After cuts: {len(data[branches[0][0]])} data points")


def read_root_branches(input_file, decay_tree, branches):
    """
    Read the provided ROOT branches into numpy arrays or something

    Returns a dict of arrays i guess

    """
    root_file = uproot4.open(input_file)
    data = root_file[decay_tree].arrays(branches, library="np")

    # Check that all our branches contain the same amount of data
    # This is the only way for our dataset to makes sense
    for key in data:
        if len(data[key]) != len(data[branches[0]]):
            raise Exception("data are not all the same length")

    return data


def main(args):
    # Locate the branches we actually might want
    # It's possible I have the 'wrong' branches here, but I'm looking for kinematic branches + branches I need to make cuts
    default_branches = [
        "Dstar_M",
        "D0_IPCHI2_OWNPV",
        "D0_PE",
        "D0_PX",
        "D0_PY",
        "D0_PZ",
        "D0_M",
        "D0_TAU",
        "Kminus_PE",
        "Kminus_PX",
        "Kminus_PY",
        "Kminus_PZ",
        "pi1plus_PE",
        "pi1plus_PX",
        "pi1plus_PY",
        "pi1plus_PZ",
        "pi2plus_PE",
        "pi2plus_PX",
        "pi2plus_PY",
        "pi2plus_PZ",
        "pi3minus_PE",
        "pi3minus_PX",
        "pi3minus_PY",
        "pi3minus_PZ",
    ]

    # If we just want to show the default branches, print and exit
    if args.show_default_branches:
        print("Default branches to read:")
        for branch in default_branches:
            print("\t", branch)
        return

    # Read the required branches into numpy arrays or something
    data = read_root_branches(
        args.inFile, args.decayTree, default_branches + args.branches
    )

    # Perform cuts
    cuts = [d_mass_cut, delta_m_cut, d_impact_param_cut]
    branches = [("D0_M",), ("Dstar_M", "D0_M"), ("D0_IPCHI2_OWNPV",)]
    perform_cuts(data, branches, cuts)

    # Read the data in to a new ROOT file
    write_root_file(data, "tmp.root", "my_tree")


def cli():
    """
    Argument parsing for this script

    Returns argparse ArgumentParser object

    """
    parser = argparse.ArgumentParser(
        description="Parsing and cutting a ROOT file down to only the required branches. Cuts are hard-coded for now"
    )
    required_arguments = parser.add_mutually_exclusive_group(required=True)
    required_arguments.add_argument("--inFile", help="ROOT file to read")
    required_arguments.add_argument(
        "--show-default-branches",
        action="store_true",
        help="Print the list of default branches and exit",
    )
    parser.add_argument(
        "--outFile", help="File to write to, defaults to 'out.root'", default="out.root"
    )
    parser.add_argument(
        "--branches",
        help="Additional branches to read. Any arguments provided are appended to the list of default branches",
        nargs="*",
        default=[],
    )
    parser.add_argument(
        "--decayTree",
        help="Decay Tree object to read, defaults to TupleDstToD0pi_D0ToKpipipi/DecayTree",
        default="TupleDstToD0pi_D0ToKpipipi/DecayTree",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main(cli())
