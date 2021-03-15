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


def test_inv_mass_consistency():
    """
    Consistency check that the invariant mass fcn using a single particle returns the same result as when using a numpy array

    """
    particle = [1, 2, 3, 4]
    assert (
        phsp_parameterisation.invariant_mass(*particle)
        == phsp_parameterisation.invariant_masses(
            np.array([particle[0]]),
            np.array([particle[1]]),
            np.array([particle[2]]),
            np.array([particle[3]]),
        )[0]
    )


def test_momentum_ordering():
    """
    Check that pions are returned in the correct order when passed to momentum_order

    """
    kaon = [100, 200, 300, 619.449]
    pi1 = [100, 200, 300, 399.349]  # Lower inv mass M(Kpi)
    pi2 = [50, -200, 150, 290.654]  # Higher inv mass M(Kpi)

    assert np.array_equal(
        phsp_parameterisation.momentum_order(kaon, pi1, pi2), (pi1, pi2)
    )
    assert np.array_equal(
        phsp_parameterisation.momentum_order(kaon, pi2, pi1), (pi1, pi2)
    )


def test_momm_ordering_array():
    """
    Check that momentum ordering works correctly when modifying numpy arrays in-place

    """
    # Create multidimensional numpy arrays of (k_px, k_py, k_pz, k_E) etc. for our particles
    # Only one particle of each species so each will be a 4x1 array
    kaon = np.array(((100,), (200,), (300,), (619.449,)))
    pi1 = np.array(((100,), (200,), (300,), (399.349,)))  # Lower inv mass M(Kpi)
    pi2 = np.array(((50,), (-200,), (150,), (290.654,)))  # Higher inv mass M(Kpi)

    # Store copies of these arrays for comparison later
    old_pi1, old_pi2 = pi1.copy(), pi2.copy()

    # Check the no-op ordering works
    pi1.T[0], pi2.T[0] = phsp_parameterisation.momentum_order(
        kaon.T[0], pi1.T[0], pi2.T[0]
    )
    assert np.array_equal(pi1, old_pi1)
    assert np.array_equal(pi2, old_pi2)

    # Check the ordering works
    pi2.T[0], pi1.T[0] = phsp_parameterisation.momentum_order(
        kaon.T[0], pi1.T[0], pi2.T[0]
    )
    assert np.array_equal(pi2, old_pi1)
    assert np.array_equal(pi1, old_pi2)
