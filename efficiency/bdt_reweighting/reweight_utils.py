"""
Stuff that's helpful for reweighting using python BDT

Python I/O, phase space parametrisation and BDT hyperparameter optimisation

"""
import uproot


def read_branch(file_name: str, tree_name: str, branch_name: str) -> np.ndarray:
    """
    Read the contents of a ROOT branch into a numpy.ndarray

    """
    tree = uproot.open(file_name)[tree_name]

    return tree.array(branch_name)


def invariant_mass(particles):
    """
    Find the invariant mass of a collection of particles

    particles should be an iterable whose data members contain .px, .py, .pz, and .energy attributes

    """
    energy = 0.0
    px, py, pz = 0.0, 0.0, 0.0

    for particle in particles:
        energy += particle.energy
        px += particle.px
        py += particle.py
        pz += particle.pz

    p_squared = px ** 2 + py ** 2 + pz ** 2
    mass_squared = energy * energy - p_squared

    return np.sqrt(mass_squared)
