import pytest
from util import util


def test_time_bin_underflow():
    bins = (1, 2, 3)

    # Check unsafe binning doesn't raise when we underflow
    util.unsafe_bin(-1, bins)

    with pytest.raises(ValueError):
        util.bin(-1.0, bins)


def test_time_bin_overflow():
    bins = (1, 2, 3)

    # Check unsafe binning doesn't raise when we overflow
    util.unsafe_bin(4, bins)

    with pytest.raises(ValueError):
        util.bin(4.0, bins)


def test_binning():
    bins = (1, 2, 3)
    assert util.bin(2.5, bins) == 1
