import pytest
from .. import categorisation


def test_time_bin_underflow():
    with pytest.raises(ValueError):
        categorisation.time_bin(-1.0)


def test_time_bin_overflow():
    with pytest.raises(ValueError):
        categorisation.time_bin(0.041)


def test_time_bin():
    time = 0.00041  # 1 D lifetime; should go in bin 2
    assert categorisation.time_bin(time) == 2
