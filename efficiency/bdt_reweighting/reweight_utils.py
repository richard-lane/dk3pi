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


def signature(k, pi1, pi2, pi3):
    """
    Returns the "signature" of a final state: i.e. sgn(det(M)), where M is a matrix [[k], [pi1], [pi2], [pi3]]

    Intended for use with kinematic data in the order (px, py, pz, E), but it doesn't really matter as long as you're consistent

    params:
      k: kaon parameters (px, py, pz, E)
      pi1: pion parameters (px, py, pz, E)
      pi2: pion parameters (px, py, pz, E)
      pi3: pion parameters (px, py, pz, E)

    returns +1 or -1

    """
    return np.sign(np.linalg.det(np.row_stack((k, pi1, pi2, pi3))))


def signatures(
    k: np.ndarray, pi1: np.ndarray, pi2: np.ndarray, pi3: np.ndarray
) -> np.ndarray:
    """
    Find the signatures of a collection of particles, where each is passed in as numpy arrays e.g.
      k = [[k_px], [k_py], [k_pz], [E]]

    """
    n = len(k[0])
    assert n == len(pi1[0]) == len(pi2[0]) == len(pi3[0])

    sigs = np.zeros(n)
    for i in range(n):
        sigs[i] = signature(k.T[i], pi1.T[i], pi2.T[i], pi3.T[i])

    return sigs


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


def read_kinematic_data(
    file_name: str, tree_name: str, k_branches, pi1_branches, pi2_branches, pi3_branches
):
    """
    Read K3pi kinematic data from a ROOT file

    branches should be iterables of branch names in the order [px, py, pz, pe]

    Returns particles K, pi1, pi2, pi3 as e.g. [[kpx...], [kpy...], [kpz...], [kE...]]

    """
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

    return k, pi1, pi2, pi3


def read_invariant_masses(
    file_name: str, tree_name: str, k_branches, pi1_branches, pi2_branches, pi3_branches
) -> np.ndarray:
    """
    Find the (Kpi1, pi1pi2, pi2pi3, Kpi1pi2, pi1pi2pi3) phsp parametrisation of a set of points in a ROOT file

    branches should be iterables of branch names in the order [px, py, pz, pe]

    Performs momentum ordering of pi1 and pi2 based on the momentum_order function

    """
    # Read the data from all the events
    k, pi1, pi2, pi3 = read_kinematic_data(
        file_name, tree_name, k_branches, pi1_branches, pi2_branches, pi3_branches
    )

    # Perform momentum ordering
    # This sometimes does assignments that it doesn't need to but it should be ok
    # TODO function this
    for i in range(len(k[0])):
        # Assign i'th pi1 and pi2 params to the right things
        pi1.T[i], pi2.T[i] = momentum_order(k.T[i], pi1.T[i], pi2.T[i])

    return invariant_mass_parametrisation(k, pi1, pi2, pi3)


def hist_error(bins, data):
    """
    Poisson errors on a binned dataset

    :param bins: The bins to consider. Should be a sorted iterable of bin edges, containing the lower edge of every bin and the higher edge of the highest one.
    :param data: The dataset to be binned

    :returns: an iterable of the Poisson errors in each bin

    """
    binned_data, _ = np.histogram(data, bins=bins)

    return np.sqrt(binned_data)


def fractional_ratio_error(bins, numerator, denominator):
    """
    Fractional errors on a ratio of histograms.

    Assumes Poisson errors for each histogram

    :param bins: The bins to consider. Should be a sorted iterable of bin edges, containing the lower edge of every bin and the higher edge of the highest one.
    :param numerator: numerator dataset
    :param denominator: numerator dataset

    :returns: an iterable of the fractional errors in each bin

    """
    numerator_binned, _ = np.histogram(numerator, bins=bins)
    denominator_binned, _ = np.histogram(denominator, bins=bins)

    # Fractional error is sqrt((n + d )/ (n * d))
    hist_sum = np.add(numerator_binned, denominator_binned)
    hist_prod = np.multiply(numerator_binned, denominator_binned)

    return np.sqrt(np.divide(hist_sum, hist_prod))
