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


def plot_projection(phsp_points, i, label):
    """
    Return a histogram of the i'th projection of a collection of phase space points

    Returns whatever plt.hist returns

    """
    # Create an array of the points we're interested in
    data = phsp_points[:, i]

    return plt.hist(data, bins=np.linspace(200, 1800, 250))


def invariant_mass_parametrisation(k, pi1, pi2, pi3):
    """
    Given arrays of
    ((k_px, k_py, k_pz, k_energy),
     (pi1_px, pi1_py, pi1_pz, pi1_energy),
     (pi2_px, pi2_py, pi2_pz, pi2_energy),
     (pi3_px, pi3_py, pi3_pz, pi3_energy))

    find the invariant mass parametrisation (kpi1, pi1pi2, pi2pi3, kpi1pi2, pi1pi2pi3)

    """
    kpi1 = invariant_masses(*np.add(k, pi1))
    pi1pi2 = invariant_masses(*np.add(pi1, pi2))
    pi2pi3 = invariant_masses(*np.add(pi2, pi3))
    kpi1pi2 = invariant_masses(*np.add(np.add(k, pi1), pi2))
    pi1pi2pi3 = invariant_masses(*np.add(np.add(pi1, pi2), pi3))

    return np.column_stack((kpi1, pi1pi2, pi2pi3, kpi1pi2, pi1pi2pi3))


def read_invariant_masses(
    file_name: str, tree_name: str, k_branches, pi1_branches, pi2_branches, pi3_branches
) -> np.ndarray:
    """
    Find the (Kpi1, pi1pi2, pi2pi3, Kpi1pi2, pi1pi2pi3) phsp parametrisation of a set of points in a ROOT file

    branches should be iterables of branch names in the order [px, py, pz, pe]

    """
    # Read the data from all the events
    k_px = read_branch(file_name, tree_name, k_branches[0])
    k_py = read_branch(file_name, tree_name, k_branches[1])
    k_pz = read_branch(file_name, tree_name, k_branches[2])
    k_e = read_branch(file_name, tree_name, k_branches[3])

    pi1_px = read_branch(file_name, tree_name, pi1_branches[0])
    pi1_py = read_branch(file_name, tree_name, pi1_branches[1])
    pi1_pz = read_branch(file_name, tree_name, pi1_branches[2])
    pi1_e = read_branch(file_name, tree_name, pi1_branches[3])

    pi2_px = read_branch(file_name, tree_name, pi2_branches[0])
    pi2_py = read_branch(file_name, tree_name, pi2_branches[1])
    pi2_pz = read_branch(file_name, tree_name, pi2_branches[2])
    pi2_e = read_branch(file_name, tree_name, pi2_branches[3])

    pi3_px = read_branch(file_name, tree_name, pi3_branches[0])
    pi3_py = read_branch(file_name, tree_name, pi3_branches[1])
    pi3_pz = read_branch(file_name, tree_name, pi3_branches[2])
    pi3_e = read_branch(file_name, tree_name, pi3_branches[3])

    return invariant_mass_parametrisation(
        (k_px, k_py, k_pz, k_e),
        (pi1_px, pi1_py, pi1_pz, pi1_e),
        (pi2_px, pi2_py, pi2_pz, pi2_e),
        (pi3_px, pi3_py, pi3_pz, pi3_e),
    )
