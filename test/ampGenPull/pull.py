import sys
import os
import subprocess
import argparse
import numpy as np


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
    parser.add_argument(
        "--out_file",
        default="pull.txt",
        type=str,
        help="file to write fit results to. format defined by the fit executable",
    )
    return parser.parse_args()


def main(args):
    # Find how many events to generate from a poisson distribution
    num_cf = args.numCF
    out_file = args.out_file
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
    dcs_root_file = "d.root"
    cf_root_file = "dBar.root"

    subprocess.run(
        [
            ampgen_generator,
            dcs_opt_file,
            "--nEvents",
            str(num_dcs),
            "--EventType",
            "D K+ pi- pi- pi+",
            "--Output",
            dcs_root_file,
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
            cf_root_file,
            "--GenerateTimeDependent",
        ],
        check=True,
    )

    # Read ROOT files, calculate ratios, fit + save the results to a file
    fitter_executable = os.path.abspath(
        os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
            "build",
            "test",
            "ampGenPull",
            "ampgenpull",
        )
    )
    subprocess.run([fitter_executable, dcs_root_file, cf_root_file], check=True)

    r_pulls = []
    with open(out_file, "r") as f:
        fit_results = [line.rstrip() for line in f]

        for result in fit_results:
            result = result.split(",")
            r = result[0]
            dr = result[1]
            r_pulls.append((float(r) - 0.0549) / float(dr))
    print(np.mean(r_pulls), "+-", np.std(r_pulls))


if __name__ == "__main__":
    main(cli())
