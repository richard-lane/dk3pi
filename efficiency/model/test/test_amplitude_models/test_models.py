import phasespace
import numpy as np

from util import definitions


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

    k = np.resize(k, num_accepted)
    pi1 = np.resize(pi1, num_accepted)
    pi2 = np.resize(pi2, num_accepted)
    pi3 = np.resize(pi3, num_accepted)

    # Weight according to CF amplitude model
    # Read AmpGen events
    # Plot projections of them both to see what they look like
    ...
