"""
Stuff that's helpful for reweighting using python BDT

Python I/O, phase space parametrisation and BDT hyperparameter optimisation

"""
import uproot
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
