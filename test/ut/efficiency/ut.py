"""
Should probably also add this to the CMake automated tests, but maybe I'll do that when I delete the C++ efficiency stuff

"""
import pytest
import numpy as np

# really
import sys
from os.path import dirname, join

sys.path.append(
    join(dirname(dirname(dirname(dirname(__file__)))), "efficiency", "bdt_reweighting")
)
import reweight_utils


def test_inv_mass():
    """
    Check that invariant mass is calculated correctly

    """
    assert reweight_utils.invariant_mass(0, 0, 0, 2) == 2
    assert np.isclose(reweight_utils.invariant_mass(100, 200, 300, 600), 469.041576)


def test_momentum_ordering():
    """
    Check that pions are returned in the correct order when passed to momentum_order

    """
    kaon = [100, 200, 300, 619.449]
    pi1 = [100, 200, 300, 399.349]  # Lower inv mass M(Kpi)
    pi2 = [50, -200, 150, 290.654]  # Higher inv mass M(Kpi)

    assert np.array_equal(reweight_utils.momentum_order(kaon, pi1, pi2), (pi1, pi2))
    assert np.array_equal(reweight_utils.momentum_order(kaon, pi2, pi1), (pi1, pi2))


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
    pi1.T[0], pi2.T[0] = reweight_utils.momentum_order(kaon.T[0], pi1.T[0], pi2.T[0])
    assert np.array_equal(pi1, old_pi1)
    assert np.array_equal(pi2, old_pi2)

    # Check the ordering works
    pi2.T[0], pi1.T[0] = reweight_utils.momentum_order(kaon.T[0], pi1.T[0], pi2.T[0])
    assert np.array_equal(pi2, old_pi1)
    assert np.array_equal(pi1, old_pi2)


def test_inv_mass_consistency():
    """
    Consistency check that the invariant mass fcn using a single particle returns the same result as when using a numpy array

    """
    particle = [1, 2, 3, 4]
    assert reweight_utils.invariant_mass(*particle) == reweight_utils.invariant_masses(
        np.array([particle[0]]),
        np.array([particle[1]]),
        np.array([particle[2]]),
        np.array([particle[3]]),
    )[0]
