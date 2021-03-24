import phasespace
import numpy as np
import matplotlib.pyplot as plt

from util import definitions
from util.phsp_parameterisation import invariant_masses
from util import phsp_binning


def mkpi1(k, pi1):
    """
    Find the masses M(K pi1) of a set of events

    """
    return invariant_masses(*np.add(k, pi1))


def test_cf():
    """
    Test CF model for consistency with AmpGen

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

    # Weight according to CF amplitude model
    definitions.assign_cf(phsp_binning._find_fcn(definitions.CF_LIB, "cf_wrapper"))
    definitions.assign_dcs(phsp_binning._find_fcn(definitions.DCS_LIB, "dcs_wrapper"))
    cf_amplitude = [
        phsp_binning._cf_amplitude(
            (np.concatenate((k.T[i], pi1.T[i], pi2.T[i], pi3.T[i]))), +1
        )
        for i in range(num_accepted)
    ]

    # Read AmpGen events

    # Plot projections of them both to see what they look like
    phsp_mkpi = mkpi1(k[:, 0:num_accepted], pi1[:, 0:num_accepted])
    plt.hist(phsp_mkpi, bins=100)
    plt.show()
