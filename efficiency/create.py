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
        times = np.concatenate((times, tree["D0_TAU"].array()))

    # Categorise events into phsp bins

    # Read AmpGen
    ag_times = np.array([])
    ag_phsp = np.array([]).reshape(0, 5)

    with uproot.open(definitions.RS_AMPGEN_PATH) as f:
        tree = f["DalitzEventList"]
        ag_times = np.concatenate((ag_times, tree["Dbar0_decayTime"].array()))

        k = [
            1000 * tree[branch].array()
            for branch in ("_1_K~_Px", "_1_K~_Py", "_1_K~_Pz", "_1_K~_E")
        ]
        pi1 = [
            1000 * tree[branch].array()
            for branch in ("_2_pi#_Px", "_2_pi#_Py", "_2_pi#_Pz", "_2_pi#_E")
        ]
        pi2 = [
            1000 * tree[branch].array()
            for branch in ("_3_pi#_Px", "_3_pi#_Py", "_3_pi#_Pz", "_3_pi#_E")
        ]
        pi3 = [
            1000 * tree[branch].array()
            for branch in ("_4_pi~_Px", "_4_pi~_Py", "_4_pi~_Pz", "_4_pi~_E")
        ]

        # Perform momentum ordering

        ag_phsp = np.concatenate(
            (
                ag_phsp,
                phsp_parameterisation.invariant_mass_parametrisation(k, pi1, pi2, pi3),
            ),
            axis=0,
        )

    # Perform Ks veto

    d_lifetime = 4.101e-4
    times /= d_lifetime
    ag_times /= d_lifetime

    times = times[times > 0]

    kw = {
        "bins": 100,
        "alpha": 0.3,
        "density": True,
    }
    plt.hist(phsp_points[:, 4], **kw)
    plt.hist(ag_phsp[:, 4], **kw)

    plt.savefig("tmp.png")

    plt.clf()
    plt.hist(times, label="MC", **kw, range=(-10, 10))
    plt.hist(ag_times, label="AmpGen", **kw, range=(-10, 10))
    plt.legend()
    plt.xlabel("D lifetimes")
    plt.savefig("times.png")


if __name__ == "__main__":
    main()
