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


def invariant_mass(px, py, pz, energy):
    """
    Find the invariant mass of a particle given its kinematic data

    """
    return np.sqrt(energy ** 2 - px ** 2 - py ** 2 - pz ** 2)


def invariant_masses(
    px: np.ndarray, py: np.ndarray, pz: np.ndarray, energy: np.ndarray
):
    """
    Find the invariant masses of a collection of particles represented by their kinematic data

    all arguments are 1d np arrays of momentum/energy

    Returns an array of invariant masses

    """
    # Turns out invariant_mass as defined above just works with numpy arrays.
    return invariant_mass(px, py, pz, energy)


def momentum_order(k, pi1, pi2):
    """
    Order two pions based on the invariant mass M(Kpi)

    params:
      k: kaon parameters (px, py, pz, E)
      pi1: pion parameters (px, py, pz, E)
      pi2: pion parameters (px, py, pz, E)

    returns (lower_mass_pion, higher_mass_pion)

    """
    pi1_copy = np.copy(pi1)
    pi2_copy = np.copy(pi2)

    m1 = invariant_mass(*np.add(k, pi1_copy))
    m2 = invariant_mass(*np.add(k, pi2_copy))
    if m1 < m2:
        return pi1_copy, pi2_copy

    return pi2_copy, pi1_copy


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

    Performs momentum ordering of pi1 and pi2 based on the momentum_order function

    """
    # Read the data from all the events
    k = np.array([read_branch(file_name, tree_name, branch) for branch in k_branches])
    pi1 = np.array(
        [read_branch(file_name, tree_name, branch) for branch in pi1_branches]
    )
    pi2 = np.array(
        [read_branch(file_name, tree_name, branch) for branch in pi2_branches]
    )
    pi3 = np.array(
        [read_branch(file_name, tree_name, branch) for branch in pi3_branches]
    )

    # Perform momentum ordering
    # This sometimes does assignments that it doesn't need to but it should be ok
    for i in range(len(k[0])):
        # Assign i'th pi1 and pi2 params to the right things
        pi1.T[i], pi2.T[i] = momentum_order(
            k.T[i],
            pi1.T[i],
            pi2.T[i],
        )

    return invariant_mass_parametrisation(k, pi1, pi2, pi3)
