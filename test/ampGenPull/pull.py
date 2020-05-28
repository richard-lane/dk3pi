import sys
import os
import subprocess
import argparse


def cli():
    """
    CLI parsing

    Returns an argparse args object
    """
    parser = argparse.ArgumentParser(
        description="Run a pull study using events generated with AmpGen"
    )
    parser.add_argument(
        "--numCF",
        default=60000,
        type=int,
        help="Number of CF events to generate. The number of DCS events is calculated using the phase space parameters used by AmpGen (defined in the relevant .opt files). Defaults to 60000",
    )
    return parser.parse_args()


def main(args):
    # Find how many events to generate from a poisson distribution
    num_cf = args.numCF
    num_dcs = int(
        num_dcs / 75
    )  # Calculated from the phase space params that we expect, given the inputs to AmpGen

    # Run AmpGen

    # Read ROOT files, calculate ratios, fit + save the results to a file
    pass


if __name__ == "__main__":
    main(cli())
