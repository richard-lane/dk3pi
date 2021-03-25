import pytest
import phasespace
from math import sqrt
import numpy as np

from util import definitions
from util import phsp_binning
from util import phsp_parameterisation


@pytest.mark.flaky
def test_models():
    """
    Test models for consistency with what we expect

    By construction, we expect phsp events to have Z = 0.47, arg(Z) = 0, r = 0.055

    """
    # Generate phsp events
    N = 100000
    max_wt = 0.12  # ish?
    generator = phasespace.nbody_decay(
        definitions.D_MASS,
        (
            definitions.K_MASS,
            definitions.PI_MASS,
            definitions.PI_MASS,
            definitions.PI_MASS,
        ),
        names=("K", "pi1", "pi2", "pi3"),
    )
    k = np.zeros((4, N))
    pi1 = np.zeros((4, N))
    pi2 = np.zeros((4, N))
    pi3 = np.zeros((4, N))

    weights, particles = generator.generate(N, normalize_weights=True)

    random_numbers = max_wt * np.random.random(N)

    num_accepted = 0
    for rnd, wt, this_k, this_pi1, this_pi2, this_pi3 in zip(
        random_numbers,
        weights,
        particles["K"].numpy(),
        particles["pi1"].numpy(),
        particles["pi2"].numpy(),
        particles["pi3"].numpy(),
    ):
        assert wt < max_wt

        if rnd < wt:
            k.T[num_accepted] = this_k
            pi1.T[num_accepted] = this_pi1
            pi2.T[num_accepted] = this_pi2
            pi3.T[num_accepted] = this_pi3

            num_accepted += 1

    num_dcs = 0.0
    num_cf = 0.0
    z = 0.0
    for this_k, this_pi1, this_pi2, this_pi3 in zip(
        k.T[0:num_accepted],
        pi1.T[0:num_accepted],
        pi2.T[0:num_accepted],
        pi3.T[0:num_accepted],
    ):
        if not phsp_parameterisation.vetoed(
            this_pi1, this_pi3
        ) and not phsp_parameterisation.vetoed(this_pi2, this_pi3):
            # Evaluate CF amplitude
            cf_amp = phsp_binning._cf_amplitude(
                (np.concatenate((this_k, this_pi1, this_pi2, this_pi3))), +1
            )

            # Evaluate DCS amplitude
            dcs_amp = (
                phsp_binning._dcs_amplitude(
                    (np.concatenate((this_k, this_pi1, this_pi2, this_pi3))), +1
                )
                * definitions.DCS_OFFSET
            )

            # Find Z for this event
            z += cf_amp * dcs_amp.conjugate()

            # Increment num_dcs and num_cf counters
            num_cf += abs(cf_amp) ** 2
            num_dcs += abs(dcs_amp) ** 2

    R = abs(z) / sqrt(num_dcs * num_cf)
    d = np.angle(z, deg=True)
    r = sqrt(num_dcs / num_cf)

    print(f"R = {R}; should be 0.47")
    print(f"d = {d}; should be 0.0")
    print(f"r = {r}; should be 0.055")

    assert abs(R - 0.47) < 0.02
    assert abs(d) < 1.0  # Less than 1 degree
    assert abs(r - 0.055) < 0.005
