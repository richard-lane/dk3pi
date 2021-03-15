import numpy as np

from util import phsp_parameterisation


def test_inv_mass():
    """
    Test invariant mass of one particle

    """
    assert phsp_parameterisation.invariant_mass(0, 0, 0, 2) == 2
    assert np.isclose(
        phsp_parameterisation.invariant_mass(100, 200, 300, 600), 469.041576
    )


def test_inv_masses():
    """
    Test invariant masses of multiple particles

    """
    k_mass = 493.677
    pi_mass = 139.570

    kaon = np.array([-294.02137063, 284.22807634, 288.87688733, 703.13654586])
    pion = np.array([-1.17766755, 288.05639323, -143.61907326, 350.83343012])
    pion2 = np.array([80.56682159, -38.96982116, 71.85998617, 180.7022474])

    particles = np.column_stack((kaon, pion, pion2))

    assert np.allclose(
        phsp_parameterisation.invariant_masses(*particles), (k_mass, pi_mass, pi_mass)
    )

