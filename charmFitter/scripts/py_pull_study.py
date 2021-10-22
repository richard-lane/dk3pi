"""
Pull study for fitter

"""
# Set up include dir
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from scipy.integrate import quad

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "build/charmFitter/python"))

import libcleoScan


def _gen_x_y(n):
    mean = (libcleoScan.WORLD_AVERAGE_X, libcleoScan.WORLD_AVERAGE_Y)
    cov = ((libcleoScan.WORLD_AVERAGE_X_ERR ** 2, libcleoScan.X_Y_CORRELATION * libcleoScan.WORLD_AVERAGE_X_ERR * libcleoScan.WORLD_AVERAGE_Y_ERR),
           (libcleoScan.X_Y_CORRELATION * libcleoScan.WORLD_AVERAGE_X_ERR * libcleoScan.WORLD_AVERAGE_Y_ERR, libcleoScan.WORLD_AVERAGE_Y_ERR ** 2))

    vals = np.random.multivariate_normal(mean, cov, size=n)
    return vals[:, 0], vals[:, 1]


def main():
    # Define our parameters
    width = 2.5
    max_time = 5.0 / width
    n_rs = 1_000_00  # Number of RS evts to generate each time
    k = 5  # Number of times to repeat generation
    num_experiments = 100
    actual_re_z, actual_r = -0.1, 0.055

    # Define bins
    bins = np.linspace(0, max_time, 10)

    # Generate X and Y values to use based on the world averages + their correlations/errors
    x_vals, y_vals = _gen_x_y(num_experiments)
    plt.hist2d(x_vals, y_vals, bins=50)
    plt.xlabel("x")
    plt.xlabel("y")
    plt.colorbar()
    plt.savefig("xyvals.png")
    plt.clf()

    # No efficiency
    efficiency = lambda x: 1

    re_z_vals, re_z_errs = [], []
    r_vals, r_errs = [], []

    for i in range(num_experiments):
         print(i, "-"*79)
         decay_params = [x_vals[i], y_vals[i], actual_r, 0.7, actual_re_z, width]

         # Create a fitter
         Fitter = libcleoScan.ConstrainedFitter(bins, decay_params, [1]*6)

         for _ in range(k):
              print("\tgenerating")
              seed = np.random.randint(0, np.iinfo(np.uint32).max)
              rs_times, ws_times = libcleoScan.simulate(n_rs, decay_params, max_time, seed)

              # Bin
              Fitter.addRSPoints(rs_times, np.ones(len(rs_times)))
              Fitter.addWSPoints(ws_times, np.ones(len(ws_times)))

         # Perform the fit
         print("fitting")
         Fitter.fixParameters(["z_im", "width"])
         result = Fitter.fit(efficiency)
         _, _, r, _, re_z, _ = result.fitParams
         _, _, r_err, _, re_z_err, _ = result.fitParamErrors

         re_z_vals.append(re_z)
         r_vals.append(r)

         re_z_errs.append(re_z_err)
         r_errs.append(r_err)

    re_z_vals = np.array(re_z_vals)
    r_vals = np.array(r_vals)

    re_z_vals = (re_z_vals - actual_re_z) / re_z_errs
    r_vals = (r_vals - actual_r) / r_errs

    print(f"{np.mean(re_z_vals)=}\t{np.std(re_z_vals)=}")
    print(f"{np.mean(r_vals)=}\t{np.std(r_vals)=}")

    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    kw = {"histtype": "step", "bins": np.linspace(-3, 3)}

    ax[0].hist(re_z_vals, **kw)
    ax[0].set_title(fr"{np.mean(re_z_vals):4.3}$\pm${np.std(re_z_vals):4.3}")
    ax[0].set_xlabel(r"$Re(Z)$")

    ax[1].hist(r_vals, **kw)
    ax[1].set_title(fr"{np.mean(r_vals):4.3}$\pm${np.std(r_vals):4.3}")
    ax[1].set_xlabel(r"$r$")
    plt.savefig("pull.png")

if __name__ == "__main__":
    main()

