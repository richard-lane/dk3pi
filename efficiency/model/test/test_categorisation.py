import pytest
from .. import categorisation


def test_time_bin_underflow():
    with pytest.raises(ValueError):
        categorisation.time_bin(-1.0)


def test_time_bin_overflow():
    with pytest.raises(ValueError):
        categorisation.time_bin(0.041)
