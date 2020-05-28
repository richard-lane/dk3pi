"""
Maybe this could do with some tests, but i probably won't bother until
something breaks

"""

import sys
import os
import subprocess
import argparse
import numpy as np


def num_events(dcs_mean, cf_mean):
    """
    Returns number of DCS and CF events to generate
    Returns (DCS counts, CF counts)

    """
    return np.random.poisson(dcs_mean), np.random.poisson(cf_mean)


def run_ampgen(options_file, root_file, num_events, event_type):
    """
    Run ampgen with the provided abspaths to options/root files

    """
    ampgen_generator = os.path.abspath(
        os.path.join(os.environ["AMPGENROOT"], "build", "bin", "Generator")
    )

    subprocess.run(
        [
            ampgen_generator,
            options_file,
            "--nEvents",
            str(num_events),
            "--EventType",
            event_type,
            "--Output",
            root_file,
            "--GenerateTimeDependent",
        ],
        check=True,
    )


def find_pull(out_file):
    """
    Find pull from data in a text file

    """
    expected_rd = 0.0549
    r_pulls = []
    with open(out_file, "r") as f:
        fit_results = [line.rstrip() for line in f]

        for result in fit_results:
            result = result.split(",")
            r = result[0]
            dr = result[1]
            r_pulls.append((float(r) - expected_rd) / float(dr))
    print(np.mean(r_pulls), "+-", np.std(r_pulls))


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
    parser.add_argument(
        "--numExperiments",
        default=1,
        type=int,
        help="number of pseudo experiments to run. defaults to 1 cus why not",
    )
    return parser.parse_args()


def main(args):
    # Calculate relative number of events from phase space params given by AmpGen
    dcs_fraction = 75

    # Filepaths and things needed for running ampgen and fit
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
    dcs_event_type = "D K+ pi- pi- pi+"
    cf_event_type = "Dbar0 K+ pi- pi- pi+"

    fitter_executable = os.path.abspath(
        os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
            "build",
            "test",
            "ampGenPull",
            "ampgenpull",
        )
    )

    # Create out_file; error if it already exists
    if os.path.exists(args.out_file):
        raise FileExistsError(args.out_file)
    open(args.out_file, "w").close()

    for _ in range(args.numExperiments):
        # Find how many events to generate from a poisson distribution
        num_dcs, num_cf = num_events(int(args.numCF / dcs_fraction), int(args.numCF))

        # Generate events
        run_ampgen(dcs_opt_file, dcs_root_file, num_dcs, dcs_event_type)
        run_ampgen(cf_opt_file, cf_root_file, num_cf, cf_event_type)

        # Read ROOT files, calculate ratios, fit + save the results to a file
        subprocess.run(
            [fitter_executable, dcs_root_file, cf_root_file, args.out_file], check=True
        )

    # Find and print pull for r
    find_pull(args.out_file)


if __name__ == "__main__":
    main(cli())
