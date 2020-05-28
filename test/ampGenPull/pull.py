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
        description="Run a pull study using events generated with AmpGen. Requires that target 'ampgenpull' is built in the 'build/test/ampGenPull' directory"
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
        num_cf / 75
    )  # Calculated from the phase space params that we expect, given the inputs to AmpGen

    # Run AmpGen
    ampgen_generator = os.path.abspath(
        os.path.join(os.environ["AMPGENROOT"], "build", "bin", "Generator")
    )
    cf_opt_file = os.path.abspath(
        os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
            "AmpGenTools",
            "options",
            "Dbar02piKpipi.opt",
        )
    )
    dcs_opt_file = os.path.abspath(
        os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
            "AmpGenTools",
            "options",
            "generate_mixing.opt",
        )
    )

    subprocess.run(
        [
            ampgen_generator,
            dcs_opt_file,
            "--nEvents",
            str(num_dcs),
            "--EventType",
            "D K+ pi- pi- pi+",
            "--Output",
            "dBar.root",
            "--GenerateTimeDependent",
        ],
        check=True,
    )
    subprocess.run(
        [
            ampgen_generator,
            cf_opt_file,
            "--nEvents",
            str(num_cf),
            "--EventType",
            "Dbar0 K+ pi- pi- pi+",
            "--Output",
            "d.root",
        ],
        check=True,
    )

    # Read ROOT files, calculate ratios, fit + save the results to a file

    pass


if __name__ == "__main__":
    main(cli())
