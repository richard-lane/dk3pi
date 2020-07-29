"""
Script for parsing our LHCb Monte Carlo dataset and just extracting the data we want/need

Implments the cuts in cuts.txt; uses uproot

"""
import uproot4


def main():
    # Locate the branches we actually might want
    # It's possible I have the 'wrong' branches here, but I'm looking for kinematic branches + branches I need to make cuts
    required_branches = {
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
    # Read these branches into a new ROOT file
    # Perform cuts
    # Write to a new file
    pass


def cli():
    """
    Argument parsing for this script

    """
    pass


if __name__ == "__main__":
    main()
