from model.creation import create_bdts
from model.creation import productions

from model.util import definitions
from model.util import phsp_parameterisation

import glob
import uproot
import numpy as np
import matplotlib.pyplot as plt


def main():
    phsp_points = np.array([]).reshape(0, 5)
    times = np.array([])

    # TODO SPLIT BY K CHARGES
    for filename in glob.glob(productions.get("MagDown", 2018, "RS")):
        # Read the file
        tree = uproot.open(filename)[definitions.RS_TREE]
        t = tree["D0_TAU"].array()

        # Find array of phsp points
        k = [
            tree[branch].array()
            for branch in ("D0_P0_PX", "D0_P0_PY", "D0_P0_PZ", "D0_P0_PE")
        ]
        pi1 = [
            tree[branch].array()
            for branch in ("D0_P1_PX", "D0_P1_PY", "D0_P1_PZ", "D0_P1_PE")
        ]
        pi2 = [
            tree[branch].array()
            for branch in ("D0_P2_PX", "D0_P2_PY", "D0_P2_PZ", "D0_P2_PE")
        ]
        pi3 = [
            tree[branch].array()
            for branch in ("D0_P3_PX", "D0_P3_PY", "D0_P3_PZ", "D0_P3_PE")
        ]

        # Perform momentum ordering

        phsp_points = np.concatenate(
            (
                phsp_points,
                phsp_parameterisation.invariant_mass_parametrisation(k, pi1, pi2, pi3),
            ),
            axis=0,
        )

        # Find array of decay times
        times = np.concatenate((times, t))

    # Categorise events into phsp bins

    # Read AmpGen
    "_1_K~_Px"
    "_1_K~_Py"
    "_1_K~_Pz"
    "_1_K~_E"
    "_2_pi#_Px"
    "_2_pi#_Py"
    "_2_pi#_Pz"
    "_2_pi#_E"
    "_3_pi#_Px"
    "_3_pi#_Py"
    "_3_pi#_Pz"
    "_3_pi#_E"
    "_4_pi~_Px"
    "_4_pi~_Py"
    "_4_pi~_Pz"
    "_4_pi~_E"

    # Perform Ks veto

    plt.hist(phsp_points[:, 0], bins=100)
    plt.savefig("tmp.png")


if __name__ == "__main__":
    main()
