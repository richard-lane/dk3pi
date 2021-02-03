"""
Should probably also add this to the CMake automated tests, but maybe I'll do that when I delete the C++ efficiency stuff

"""
import pytest
import numpy as np
import warnings

# really
import sys
from os.path import dirname, join

sys.path.append(
    join(dirname(dirname(dirname(dirname(__file__)))), "efficiency", "bdt_reweighting")
)
sys.path.append(
    join(dirname(dirname(dirname(dirname(__file__)))), "efficiency", "scripts")
)
import reweight_utils
import script_util


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
    assert (
        reweight_utils.invariant_mass(*particle)
        == reweight_utils.invariant_masses(
            np.array([particle[0]]),
            np.array([particle[1]]),
            np.array([particle[2]]),
            np.array([particle[3]]),
        )[0]
    )


def test_signature():
    """
    Check that the signature of a decay is correctly identified

    """
    assert (
        reweight_utils.signature([1, 3, 5, 9], [1, 3, 1, 7], [4, 3, 9, 7], [5, 2, 0, 9])
        == -1
    )
    assert (
        reweight_utils.signature([9, 3, 5, 1], [7, 3, 1, 1], [7, 3, 9, 4], [9, 2, 0, 5])
        == 1
    )


def test_signatures():
    """
    Check that the signatures of a collection of decays are correctly identified

    """
    k = np.array([[1, 9], [3, 3], [5, 5], [9, 1]])
    pi1 = np.array([[1, 7], [3, 3], [1, 1], [7, 1]])
    pi2 = np.array([[4, 7], [3, 3], [9, 9], [7, 4]])
    pi3 = np.array([[5, 9], [2, 2], [0, 0], [9, 5]])

    assert reweight_utils.signatures(k, pi1, pi2, pi3)[0] == -1
    assert reweight_utils.signatures(k, pi1, pi2, pi3)[1] == 1


def test_unweighted_chisq():
    """
    Check that an unweighted histogram returns the correct chi squared

    Chisq in this case should be Sum((xi-yi)^2/(xi + yi))

    """
    bins = [1, 2, 3, 4]  # 3 bins

    x = [1.5, 1.5, 2.5, 3.5, 3.5, 3.5]  # Bin content (2, 1, 3)
    y = [1.5, 2.5, 2.5, 2.5, 3.5, 3.5]  # Bin content (1, 3, 2)

    expected_chisq = 23.0 / 15.0  # Expect chisq 1.533333
    expected_p = 0.6746

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        chisq, p = script_util.chi_sq_distance(x, y, bins)

    assert np.isclose(chisq, expected_chisq, atol=1e-6)
    assert np.isclose(p, expected_p, atol=1e-4)


def test_weighted_chisq():
    """
    Check that we correctly identify the chisq distance between two weighted histograms

    """
    bins = [1, 2, 3, 4]  # 3 bins

    x = [1.5, 1.5, 2.5, 3.5, 3.5, 3.5]
    w_x = [1, 0.5, 0.5, 1, 0.5, 1]  # Bin content (1.5, 0.5, 2.5)

    y = [1.5, 2.5, 2.5, 2.5, 3.5, 3.5]
    w_y = [1, 0.5, 0.5, 1, 0.5, 1]  # Bin content (1, 2, 1.5)

    expected_chisq = 106.0 / 63.0  # Expect chisq 1.682539683
    expected_p = 0.6408

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        chisq, p = script_util.chi_sq_distance(x, y, bins, x_weights=w_x, y_weights=w_y)

    assert np.isclose(chisq, expected_chisq, atol=1e-6)
    assert np.isclose(p, expected_p, atol=1e-4)


def test_chisq_wrong_dimension():
    """
    Check that a ValueError is raised if the shape of weights and points doesn't match

    """
    bins = [1, 2, 3, 4]  # 3 bins

    good = [1.5, 1.5, 2.5, 3.5, 3.5, 3.5]
    good_weight = [1, 0.5, 0.5, 1, 0.5, 1]

    bad = [1.5, 2.5, 2.5, 2.5, 3.5, 3.5]
    bad_weight = [1, 0.5, 0.5, 1, 0.5]

    with pytest.raises(ValueError):
        script_util.chi_sq_distance(good, bad, bins, good_weight, bad_weight)

    with pytest.raises(ValueError):
        script_util.chi_sq_distance(bad, good, bins, bad_weight, good_weight)

    with pytest.raises(ValueError):
        script_util.chi_sq_distance(bad, bad, bins, bad_weight, bad_weight)


def test_chisq_underflow():
    """
    Check that a ValueError is raised if a point underflows

    """
    bins = [1, 2, 3, 4]

    good = [1.5, 1.5, 2.5, 3.5, 3.5, 3.5]

    bad = [0.5, 2.5, 2.5, 2.5, 3.5, 3.5]

    with pytest.raises(ValueError):
        script_util.chi_sq_distance(good, bad, bins)

    with pytest.raises(ValueError):
        script_util.chi_sq_distance(bad, good, bins)

    with pytest.raises(ValueError):
        script_util.chi_sq_distance(bad, bad, bins)


def test_chisq_overflow():
    """
    Check that a ValueError is raised if a point overflows

    """
    bins = [1, 2, 3, 4]

    good = [1.5, 1.5, 2.5, 3.5, 3.5, 3.5]

    bad = [2.5, 2.5, 2.5, 3.5, 3.5, 4.5]

    with pytest.raises(ValueError):
        script_util.chi_sq_distance(good, bad, bins)

    with pytest.raises(ValueError):
        script_util.chi_sq_distance(bad, good, bins)

    with pytest.raises(ValueError):
        script_util.chi_sq_distance(bad, bad, bins)

