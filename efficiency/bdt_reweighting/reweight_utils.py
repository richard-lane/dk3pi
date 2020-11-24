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

    # Find its invariant masses
    k_masses = invariant_masses(k_px, k_py, k_pz, k_e)
    phsp = np.column_stack((k_masses, k_masses))

    # Plot projections
    plot_projection(phsp, 0)


if __name__ == "__main__":
    main()
