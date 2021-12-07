"""
Plot CLEO likelihoods using Python

"""

# Set up include dir
import numpy as np
import sys, os
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), "build/charmFitter/python"))

import cleoScan

def main():
    # Create arrays from -1 to +1 for Re and Im Z
    n_z_vals = 100
    re_z = np.linspace(-1, 1, n_z_vals)
    im_z = np.linspace(-1, 1, n_z_vals)

    # Assign nonsense values of likelihood to 0.0
    nonsense = 0.0

    decay_params = [0.03, 0.06, 0.055, 0.5, 0.5, 2.5]
    tmp = cleoScan.cleoLikelihoods(re_z, im_z, np.array(decay_params), 3)
    tmp *= -2
    tmp[np.isnan(tmp)] = nonsense

    l = np.zeros((n_z_vals, n_z_vals))
    for r in range(n_z_vals):
        for i in range(n_z_vals):
            l[i, r] = tmp[i + r * n_z_vals]

    # Find min likelihood
    min_l = np.min(l[np.nonzero(l)])

    # Remove the min likelihood from each value
    l[np.nonzero(l)] -= min_l

    # Take the sqrt to transform into std dev
    l[np.nonzero(l)] = np.sqrt(l[np.nonzero(l)])

    # Plot histogram of likelihoods
    plt.hist(l.flatten(), bins=100)
    plt.savefig("py_hist.png")
    plt.clf()

    # Contour plot of likelihoods
    plt.set_cmap("jet")
    plt.pcolormesh(re_z, im_z, l)
    plt.colorbar()
    plt.savefig("py_cleo.png")


if __name__ == "__main__":
    main()

