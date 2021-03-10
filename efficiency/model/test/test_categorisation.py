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


def test_veto():
    accepted = 0.6, 0.2
    rejected = 0.49, 0.498

    for m in accepted:
        assert not categorisation.vetoed(m)

    for m in rejected:
        assert categorisation.vetoed(m)
