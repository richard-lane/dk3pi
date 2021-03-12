"""
The phase space parameterisation used to find BDT training variables

"""
from math import sqrt
import numpy as np


def invariant_mass(px, py, pz, energy):
    """
    Find the invariant mass of a particle given its kinematic data
    Units depend on input units; will be e.g. mass in MeV if momentum/energy entered in MeV

    :param px: particle x momentum
    :param py: particle y momentum
    :param pz: particle z momentum
    :param energy: particle energy
    :returns: invariant mass of the particle

    """
    # Math sqrt slightly faster than np sqrt for
    return sqrt(energy ** 2 - px ** 2 - py ** 2 - pz ** 2)


def invariant_masses(
    px: np.ndarray, py: np.ndarray, pz: np.ndarray, energy: np.ndarray
):
    """
    Find the invariant masses of a collection of particles represented by their kinematic data

    :param px: particle x momenta
    :param py: particle y momenta
    :param pz: particle z momenta
    :param energy: particle energies
    :returns: array of particle invariant masses

    """
    return np.sqrt(energy ** 2 - px ** 2 - py ** 2 - pz ** 2)


def momentum_order(k, pi1, pi2):
    """
    Order two pions based on the invariant mass M(Kpi)

    :param k: kaon parameters (px, py, pz, E)
    :param pi1: pion parameters (px, py, pz, E)
    :param pi2: pion parameters (px, py, pz, E)

    :returns: (lower_mass_pion, higher_mass_pion) as their pion parameters. Returns copies of the original arguments

    """
    pi1_copy = np.copy(pi1)
    pi2_copy = np.copy(pi2)

    m1 = invariant_mass(*np.add(k, pi1_copy))
    m2 = invariant_mass(*np.add(k, pi2_copy))
    if m1 < m2:
        return pi1_copy, pi2_copy

    return pi2_copy, pi1_copy


def invariant_mass_parametrisation(k, pi1, pi2, pi3):
    """
    Phase space parametrisation (kpi1, pi1pi2, pi2pi3, kpi1pi2, pi1pi2pi3)

    :param k: array of k parameters (k_px, k_py, k_pz, k_energy)
    :param pi1: array of pion parameters (pi1_px, pi1_py, pi1_pz, pi1_energy)
    :param pi2: array of pion parameters (pi2_px, pi2_py, pi2_pz, pi2_energy)
    :param pi2: array of pion parameters (pi3_px, pi3_py, pi3_pz, pi3_energy)

    :returns: numpy array of phase space points according to the above parameterisation.
              Shape N, 5

    """
    kpi1 = invariant_masses(*np.add(k, pi1))
    pi1pi2 = invariant_masses(*np.add(pi1, pi2))
    pi2pi3 = invariant_masses(*np.add(pi2, pi3))
    kpi1pi2 = invariant_masses(*np.add(np.add(k, pi1), pi2))
    pi1pi2pi3 = invariant_masses(*np.add(np.add(pi1, pi2), pi3))

    return np.column_stack((kpi1, pi1pi2, pi2pi3, kpi1pi2, pi1pi2pi3))
