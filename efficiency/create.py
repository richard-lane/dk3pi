from model.creation import create_bdts
from model.creation import productions

from model.util import definitions
from model.util import phsp_parameterisation
from model.util import phsp_binning

import glob
import uproot
import numpy as np
import matplotlib.pyplot as plt


def _plots(mc_points, mc_times, ag_points, ag_times):
    """
    Make some plots to test things

    """
    d_lifetime = 4.101e-4
    mc_times /= d_lifetime
    ag_times /= d_lifetime

    _, ax = plt.subplots(2, 3)

    kw = {
        "bins": 100,
        "alpha": 0.3,
        "density": True,
    }

    ax[0][0].hist(mc_times, label="mc", **kw, range=(-10, 10))
    ax[0][0].hist(ag_times, label="AmpGen", **kw, range=(-10, 10))
    ax[0][0].set_title("Decay Time")

    for i, j, t in (
        (0, 1, r"$M(K\pi_1)$"),
        (0, 2, r"$M(\pi_1\pi_2)$"),
        (1, 0, r"$M(\pi_2\pi_3)$"),
        (1, 1, r"$M(K\pi_1\pi_2)$"),
        (1, 2, r"$M(\pi_1\pi_2\pi_3)$"),
    ):
        ax[i][j].hist(mc_points[:, i], **kw)
        ax[i][j].hist(ag_points[:, i], **kw)
        ax[i][j].set_title(t)

    plt.savefig("tmp.png")


def _phsp_param(tree, kinematic_branches, time_branch, units):
    """
    Read the kinematic branches from the provided TTree and use them to find the phsp parametrisation.

    Finds the phsp parametrisation of the events in the tree, removing events close to the Ks resonance according
    to the veto criteria in util.definitions.
    Orders the same-charge pions pi1, pi2 according to util.phsp_parameterisation.momentum_ordering

    :param tree: TTree to read data from
    :param kinematic_branches: 16 k, pi1, pi2, pi2 branch names in the order (k_px, kpy, kpz, ke...) etc for k, pi1, pi2, pi3
    :param time_branch: branch name for the decay time
    :param units: units that data is in. Should be one of MeV/GeV

    :returns: phsp parametrisation according to util.phsp_parametrisation.inv_mass_parametrisaion, in GeV
    :returns: decay times
    :returns: array of phsp bin numbers, using the amplitude models in util.phsp_binning

    """
    assert units in {"MeV", "GeV"}
    assert len(kinematic_branches) == 16

    # Read branches
    k = np.array([tree[branch].array() for branch in kinematic_branches[0:4]])
    pi1 = np.array([tree[branch].array() for branch in kinematic_branches[4:8]])
    pi2 = np.array([tree[branch].array() for branch in kinematic_branches[8:12]])
    pi3 = np.array([tree[branch].array() for branch in kinematic_branches[12:16]])

    # We might have to convert to GeV
    if units == "MeV":
        k /= 1000.0
        pi1 /= 1000.0
        pi2 /= 1000.0
        pi3 /= 1000.0

    # Perform momentum ordering
    for i in range(len(k.T)):
        pi1.T[i], pi2.T[i] = phsp_parameterisation.momentum_order(
            k.T[i], pi1.T[i], pi2.T[i]
        )

    # Find phsp parameterisation
    phsp_points = phsp_parameterisation.invariant_mass_parametrisation(k, pi1, pi2, pi3)

    # Find decay time
    decay_times = tree[time_branch].array()

    # Find phsp bin for each point
    events = [
        np.concatenate((k.T[i], pi1.T[i], pi2.T[i], pi3.T[i])) for i in range(len(k.T))
    ]
    bins = np.array([phsp_binning.phsp_bin(event, +1) for event in events])

    # Find which indices need to be thrown away because of the Ks veto
    pi1pi3_masses = phsp_parameterisation.invariant_masses(*(np.add(pi1, pi3)))
    pi2pi3_masses = phsp_parameterisation.invariant_masses(*(np.add(pi2, pi3)))

    delete_indices = np.where(
        np.logical_and(
            abs(pi1pi3_masses - definitions.KS_MASS) < definitions.VETO_WIDTH,
            abs(pi2pi3_masses - definitions.KS_MASS) < definitions.VETO_WIDTH,
        )
    )

    # Throw away the points we don't need
    return (
        np.delete(phsp_points, delete_indices, axis=0),
        np.delete(decay_times, delete_indices),
        np.delete(bins, delete_indices),
    )


def main():
    phsp_points = np.array([]).reshape(0, 5)
    times = np.array([])
    phsp_bin = np.array([])  # Which phsp bin each phsp point belongs in

    # TODO SPLIT BY K CHARGES
    for filename in glob.glob(productions.get("MagDown", 2018, "RS")):
        tree = uproot.open(filename)[definitions.RS_TREE]
        phsp, t, bins = _phsp_param(
            tree,
            (
                "D0_P0_PX",
                "D0_P0_PY",
                "D0_P0_PZ",
                "D0_P0_PE",
                "D0_P1_PX",
                "D0_P1_PY",
                "D0_P1_PZ",
                "D0_P1_PE",
                "D0_P2_PX",
                "D0_P2_PY",
                "D0_P2_PZ",
                "D0_P2_PE",
                "D0_P3_PX",
                "D0_P3_PY",
                "D0_P3_PZ",
                "D0_P3_PE",
            ),
            "D0_TAU",
            "MeV",
        )

        phsp_points = np.concatenate((phsp_points, phsp), axis=0)
        times = np.concatenate((times, t))
        phsp_bin = np.concatenate((phsp_bin, bins))

    # Read AmpGen
    ampgen_tree = uproot.open(
        "/afs/cern.ch/user/j/jsmallwo/public/ForRichard/AmpGen/generator_mixing_RS_TD_100000.root"
    )["DalitzEventList"]
    ag_phsp_points, ag_times, ag_bins = _phsp_param(
        ampgen_tree,
        (
            "_1_K~_Px",
            "_1_K~_Py",
            "_1_K~_Pz",
            "_1_K~_E",
            "_2_pi#_Px",
            "_2_pi#_Py",
            "_2_pi#_Pz",
            "_2_pi#_E",
            "_3_pi#_Px",
            "_3_pi#_Py",
            "_3_pi#_Pz",
            "_3_pi#_E",
            "_4_pi~_Px",
            "_4_pi~_Py",
            "_4_pi~_Pz",
            "_4_pi~_E",
        ),
        "D_decayTime",
        "GeV",
    )

    _plots(phsp_points, times, ag_phsp_points, ag_times)


if __name__ == "__main__":
    main()
