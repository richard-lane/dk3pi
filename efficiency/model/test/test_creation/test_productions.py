from ...creation import productions
import pytest


def test_get_production():
    """
    Test that we can uniquely identify an existing analysis production using decay type, year and magnetisation

    """
    assert productions.get("MagDown", 2018, "RS")


def test_get_nonexistent_production():
    """
    Test that attempting to get a nonexistent production throws a KeyError

    """
    with pytest.raises(KeyError):
        productions.get("MagDown", 2015, "RS")


def test_get_production_bad_args():
    """
    Test that attempting to get a production with bad args throws an AssertionError

    """
    with pytest.raises(AssertionError):
        productions.get("MagSideways", 2018, "RS")

    with pytest.raises(AssertionError):
        productions.get("MagDown", "Tuesday", "RS")

    with pytest.raises(AssertionError):
        productions.get("MagDown", 2018, "no")
