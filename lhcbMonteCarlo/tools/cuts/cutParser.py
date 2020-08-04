"""
Script for reading a module of cuts implemented in python and applying them to a dataset

Needs to import a module from the command line, find what cuts it contains and apply them to our dataset

"""
import argparse
import importlib


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
        "--decayTree",
        help="Decay Tree object to read, defaults to TupleDstToD0pi_D0ToKpipipi/DecayTree for no reason",
        default="TupleDstToD0pi_D0ToKpipipi/DecayTree",
    )

    return parser.parse_args()


def main(args):
    cuts_lib = importlib.import_module(args.cuts_module)
    cuts_lib.not_a_cut_function(2)


if __name__ == "__main__":
    main(cli())
