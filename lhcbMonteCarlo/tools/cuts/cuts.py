"""
Script for parsing our LHCb Monte Carlo dataset and just extracting the data we want/need

Implments the cuts in cuts.txt; uses uproot

"""
import argparse
import uproot4


def main(args):
    # Locate the branches we actually might want
    # It's possible I have the 'wrong' branches here, but I'm looking for kinematic branches + branches I need to make cuts
    # Read these branches into a new ROOT file
    # Perform cuts
    # Write to a new file
    default_branches = {
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
        "pi2Plus_PY",
        "pi2Plus_PZ",
        "pi3minus_PE",
        "pi3minus_PX",
        "pi3minus_PY",
        "pi3minus_PZ",
    }

    # If we just want to show the default branches, then print and exit
    if args.show_default_branches:
        print("Default branches to read:")
        for branch in default_branches:
            print("\t", branch)
        return

    # Set of all our branches
    branches = default_branches.union(args.branches)


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

    parser.add_argument("--outFile", help="File to write to, defaults to 'out.root'", default="out.root")

    parser.add_argument(
        "--branches",
        help="Additional branches to read. Any arguments provided are appended to the list of default branches",
        nargs="*",
        default={},
    )

    return parser.parse_args()


if __name__ == "__main__":
    main(cli())
