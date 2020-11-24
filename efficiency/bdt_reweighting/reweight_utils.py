"""
Stuff that's helpful for reweighting using python BDT

Python I/O, phase space parametrisation and BDT hyperparameter optimisation

"""
import uproot
import numpy as np
import matplotlib.pyplot as plt


def read_branch(file_name: str, tree_name: str, branch_name: str) -> np.ndarray:
    """
    Read the contents of a ROOT branch into a numpy.ndarray

    """
    tree = uproot.open(file_name)[tree_name]

    return tree.array(branch_name)


def invariant_masses(px, py, pz, energy):
    """
    Find the invariant masses of a collection of particles represented by their kinematic data

    all arguments are 1d iterables of momentum/energy

    """
    assert len(px) == len(py) == len(pz) == len(energy)

    p_squared = px ** 2 + py ** 2 + pz ** 2
    mass_squared = energy ** 2 - p_squared

    return np.sqrt(mass_squared)


def plot_projection(phsp_points, i):
    """
    Plot a histogram of the i'th projection of a collection of phase space points

    """
    # Create an array of the points we're interested in
    data = phsp_points[:, i]

    plt.hist(data, bins=100)
    plt.show()


def main():
    """
    Test

    """
    # Read in kinematic data from a ROOT file
    f = "wg_rs_prompt.root"
    k_px = read_branch(f, "DecayTree", "D0_P0_PX")
    k_py = read_branch(f, "DecayTree", "D0_P0_PY")
    k_pz = read_branch(f, "DecayTree", "D0_P0_PZ")
    k_e = read_branch(f, "DecayTree", "D0_P0_PE")

    pi1_px = read_branch(f, "DecayTree", "D0_P1_PX")
    pi1_py = read_branch(f, "DecayTree", "D0_P1_PY")
    pi1_pz = read_branch(f, "DecayTree", "D0_P1_PZ")
    pi1_e = read_branch(f, "DecayTree", "D0_P1_PE")

    pi2_px = read_branch(f, "DecayTree", "D0_P2_PX")
    pi2_py = read_branch(f, "DecayTree", "D0_P2_PY")
    pi2_pz = read_branch(f, "DecayTree", "D0_P2_PZ")
    pi2_e = read_branch(f, "DecayTree", "D0_P2_PE")

    pi3_px = read_branch(f, "DecayTree", "D0_P3_PX")
    pi3_py = read_branch(f, "DecayTree", "D0_P3_PY")
    pi3_pz = read_branch(f, "DecayTree", "D0_P3_PZ")
    pi3_e = read_branch(f, "DecayTree", "D0_P3_PE")

    # Find its invariant masses
    k_pi1 = invariant_masses(k_px + pi1_px, k_py + pi1_py, k_pz + pi1_pz, k_e + pi1_e)
    pi1_pi2 = invariant_masses(
        pi1_px + pi2_px, pi1_py + pi2_py, pi1_pz + pi2_pz, pi1_e + pi2_e
    )
    pi2_pi3 = invariant_masses(
        pi2_px + pi3_px, pi2_py + pi3_py, pi2_pz + pi3_pz, pi2_e + pi3_e
    )
    k_pi1_pi2 = invariant_masses(
        k_px + pi1_px + pi2_px,
        k_py + pi1_py + pi2_py,
        k_pz + pi1_pz + pi2_pz,
        k_e + pi1_e + pi2_e,
    )
    pi1_pi2_pi3 = invariant_masses(
        pi1_px + pi2_px + pi3_px,
        pi1_py + pi2_py + pi3_py,
        pi1_pz + pi2_pz + pi3_pz,
        pi1_e + pi2_e + pi3_e,
    )

    points = np.column_stack((k_pi1, pi1_pi2, pi2_pi3, k_pi1_pi2, pi1_pi2_pi3))

    # Plot projections
    for i in range(5):
        plot_projection(points, i)


if __name__ == "__main__":
    main()
