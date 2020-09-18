"""
Test BDT efficiency using the Python API directly

"""

import pytest
import numpy as np
from flaky import flaky

import reweighting


@flaky
def test_bdt_flat_efficiency():
    """
    Test that our BDT can correctly identify a flat efficiency function

    """
    # Generate loads of events
    # Let's say our events are normally distributed in each dimension with a mean of 5 and std dev of 1
    # Test in 1d for Speed
    num_events = 10 ** 5
    dimensionality = 1
    mean = 5
    std_dev = 1
    generated_events = np.random.normal(
        mean, std_dev, size=(num_events, dimensionality)
    )
    detected_events = np.random.normal(
        mean, std_dev, size=(int(num_events / 2), dimensionality)
    )

    # Train BDT
    bdt = reweighting.init(generated_events, detected_events)

    # Check that the efficiency is correct
    predicted_weight = reweighting.predicted_weights(bdt, detected_events, num_events)
    assert pytest.approx(predicted_weight[0], 0.1) == 2
