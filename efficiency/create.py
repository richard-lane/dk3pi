from model.creation import create_bdts
from model.creation import productions

from model.util import definitions
from model.util import phsp_parameterisation
from model.util import phsp_binning

import glob
import uproot
import numpy as np
import matplotlib.pyplot as plt


def main():
    phsp_points = np.array([]).reshape(0, 5)
    phsp_bin = np.array([])  # Which phsp bin each phsp point belongs in
    times = np.array([])

    # TODO SPLIT BY K CHARGES
    for filename in glob.glob(productions.get("MagDown", 2018, "RS")):
        # Read the file
        tree = uproot.open(filename)[definitions.RS_TREE]

        # Find array of phsp points
        k = np.array(
            [
                0.001 * tree[branch].array()
                for branch in ("D0_P0_PX", "D0_P0_PY", "D0_P0_PZ", "D0_P0_PE")
            ]
        )
        pi1 = np.array(
            [
                0.001 * tree[branch].array()
                for branch in ("D0_P1_PX", "D0_P1_PY", "D0_P1_PZ", "D0_P1_PE")
            ]
        )
        pi2 = np.array(
            [
                0.001 * tree[branch].array()
                for branch in ("D0_P2_PX", "D0_P2_PY", "D0_P2_PZ", "D0_P2_PE")
            ]
        )
        pi3 = np.array(
            [
                0.001 * tree[branch].array()
                for branch in ("D0_P3_PX", "D0_P3_PY", "D0_P3_PZ", "D0_P3_PE")
            ]
        )

        # Perform momentum ordering
        # Should probably do this with the native arrays TODO
        for i in range(len(k)):
            pi1.T[i], pi2.T[i] = phsp_parameterisation.momentum_order(
                k.T[i], pi1.T[i], pi2.T[i]
            )

        phsp_points = np.concatenate(
            (
                phsp_points,
                phsp_parameterisation.invariant_mass_parametrisation(k, pi1, pi2, pi3),
            ),
            axis=0,
        )

        pi1pi2_masses = phsp_parameterisation.invariant_masses(*(np.add(pi1, pi2)))
        pi2pi3_masses = phsp_parameterisation.invariant_masses(*(np.add(pi2, pi3)))
        mask = np.logical_and(
            abs(pi1pi2_masses - definitions.KS_MASS) < definitions.VETO_WIDTH,
            abs(pi2pi3_masses - definitions.KS_MASS) < definitions.VETO_WIDTH,
        )

        phsp_points = np.delete(phsp_points, np.where(mask), axis=0)

        # Find array of decay times
        times = np.concatenate((times, tree["D0_TAU"].array()))
        times = np.delete(phsp_points, np.where(mask))

        # Find which phsp bin each event belongs in
        events = [
            np.concatenate((k.T[i], pi1.T[i], pi2.T[i], pi3.T[i]))
            for i in range(len(k.T))
        ]
        bins = np.array([phsp_binning.phsp_bin(event, +1) for event in events])
        phsp_bin = np.concatenate((phsp_bin, bins))
        phsp_bin = np.delete(phsp_bin, np.where(mask))

    print(phsp_bin)

    # Categorise events into phsp bins

    # Read AmpGen
    ag_times = None
    ag_phsp = None
    ag_bins = None

    with uproot.open(definitions.RS_AMPGEN_PATH) as f:
        tree = f["DalitzEventList"]

        k = np.array(
            [
                tree[branch].array()
                for branch in ("_1_K~_Px", "_1_K~_Py", "_1_K~_Pz", "_1_K~_E")
            ]
        )
        pi1 = np.array(
            [
                tree[branch].array()
                for branch in ("_2_pi#_Px", "_2_pi#_Py", "_2_pi#_Pz", "_2_pi#_E")
            ]
        )
        pi2 = np.array(
            [
                tree[branch].array()
                for branch in ("_3_pi#_Px", "_3_pi#_Py", "_3_pi#_Pz", "_3_pi#_E")
            ]
        )
        pi3 = np.array(
            [
                tree[branch].array()
                for branch in ("_4_pi~_Px", "_4_pi~_Py", "_4_pi~_Pz", "_4_pi~_E")
            ]
        )

        for i in range(len(k)):
            pi1.T[i], pi2.T[i] = phsp_parameterisation.momentum_order(
                k.T[i], pi1.T[i], pi2.T[i]
            )

        ag_phsp = phsp_parameterisation.invariant_mass_parametrisation(k, pi1, pi2, pi3)

        # Perform Ks veto
        # Could speed this up
        pi1pi2_masses = phsp_parameterisation.invariant_masses(*(np.add(pi1, pi2)))
        pi2pi3_masses = phsp_parameterisation.invariant_masses(*(np.add(pi2, pi3)))
        mask = np.logical_and(
            abs(pi1pi2_masses - definitions.KS_MASS) < definitions.VETO_WIDTH,
            abs(pi2pi3_masses - definitions.KS_MASS) < definitions.VETO_WIDTH,
        )
        ag_phsp = np.delete(ag_phsp, np.where(mask), axis=0)

        ag_times = tree["Dbar0_decayTime"].array()
        ag_times = np.delete(ag_times, np.where(mask))

        events = [
            np.concatenate((k.T[i], pi1.T[i], pi2.T[i], pi3.T[i]))
            for i in range(len(k.T))
        ]
        ag_bins = np.array([phsp_binning.phsp_bin(event, +1) for event in events])
        ag_bins = np.delete(phsp_bin, np.where(mask))

    d_lifetime = 4.101e-4
    times /= d_lifetime
    ag_times /= d_lifetime

    times = times[times > 0]

    kw = {
        "bins": 100,
        "alpha": 0.3,
        "density": False,
    }
    plt.hist(phsp_points[:, 2], **kw)
    plt.hist(ag_phsp[:, 2], **kw)
    plt.hist(
        np.linspace(
            definitions.KS_MASS - definitions.VETO_WIDTH,
            definitions.KS_MASS + definitions.VETO_WIDTH,
            100000,
        ),
        **kw
    )

    plt.savefig("tmp.png")

    plt.clf()
    plt.hist(times, label="MC", **kw, range=(-10, 10))
    plt.hist(ag_times, label="AmpGen", **kw, range=(-10, 10))
    plt.legend()
    plt.xlabel("D lifetimes")
    plt.savefig("times.png")


if __name__ == "__main__":
    main()
