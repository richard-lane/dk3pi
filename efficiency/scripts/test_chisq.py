"""
Manual test of chi squared distance

I dont think these tests are implemented right

"""
import sys
import os
import numpy as np
import phasespace
import matplotlib.pyplot as plt
from scipy.stats import chi2
from tqdm import tqdm

sys.path.append(
    os.path.abspath(
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "bdt_reweighting")
    )
)
import script_util


def _find_chisqs(num_experiments, bins, generator):
    """
    Find the chi squared distance + p values for two 1d histograms

    generator should be a function that returns two histograms and a set of weights for the first hist

    """
    chisqs = np.zeros(num_experiments)
    p_vals = np.zeros(num_experiments)

    with tqdm(total=num_experiments) as pbar:
        for i in range(num_experiments):
            a, b, weights = generator()
            chisqs[i], p_vals[i] = script_util.chi_sq_distance(
                a, b, bins, x_weights=weights
            )
            pbar.update(1)

    return chisqs, p_vals


def _plot(a, b, a_weight, bins, chisqs, p_vals):
    _, ax = plt.subplots(1, 3)

    # Plot an example of the histograms
    kw = {"bins": bins, "alpha": 0.3, "density": True}
    ax[0].hist(a, weights=a_weight, **kw)
    ax[0].hist(b, **kw)
    ax[0].set_title("Example hist")

    # Plot the distribution of chi squareds
    num_bins = len(bins) - 1
    chisq_range = np.linspace(0, 3 * num_bins)
    expected_chisq = [chi2.pdf(x, num_bins) for x in chisq_range]
    ax[1].hist(chisqs, bins=chisq_range, density=True)
    ax[1].plot(chisq_range, expected_chisq)
    ax[1].set_title("Chisq for test distribution")

    # Plot the cumulative p value
    p_vals.sort()
    ax[2].step(p_vals, np.linspace(0, num_bins, num=p_vals.size))
    ax[2].grid()
    ax[2].set_title("Cumulative p Values")
    plt.show()


def test_simple_chisq():
    """
    Check that two samples from the same distribution look the same

    """
    num_points = 2500
    num_experiments = 100
    num_bins = 50
    bins = np.linspace(0, 1, num=num_bins + 1)

    def gen():
        a, b = np.split(np.random.rand(2 * num_points), 2)
        return a, b, None

    chisqs, p_vals = _find_chisqs(num_experiments, bins, gen)

    _plot(*gen(), bins, chisqs, p_vals)


def test_weighted_chisq():
    """
    Check that we get the right chisq distribution when taking chisq between a weighted and unweighted distribution

    """
    num_points = 1000
    num_experiments = 100
    num_bins = 20
    bins = np.linspace(-800, 800, num=num_bins + 1)

    pi_mass = 139.570
    k_mass = 493.677
    d_mass = 1864.84
    generator = phasespace.nbody_decay(
        d_mass, (k_mass, pi_mass, pi_mass, pi_mass), names=("K", "pi1", "pi2", "pi3")
    )

    def gen():
        """
        Return K_px arrays a and b for a D->K3pi event, and weights for a

        Returns a, b, wt_a

        """
        # Find k_px with weights
        a_wt, a = generator.generate(num_points)
        a = a["K"].numpy()[:, 0]

        # Normalise weights to have an average of 1
        a_wt /= np.mean(a_wt)

        # Find k_px using accept-reject
        b = script_util.flat_phsp_points(num_points)[0][0]

        return a, b, a_wt

    chisqs, p_vals = _find_chisqs(num_experiments, bins, gen)
    _plot(*gen(), bins, chisqs, p_vals)


if __name__ == "__main__":
    test_simple_chisq()
    test_weighted_chisq()
