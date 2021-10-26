"""
Pull study for fitter

"""
# Set up include dir
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from multiprocessing import Queue, Process

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "build/charmFitter/python"))

import libcleoScan


def _gen_x_y(n):
    mean = (libcleoScan.WORLD_AVERAGE_X, libcleoScan.WORLD_AVERAGE_Y)
    cov = ((libcleoScan.WORLD_AVERAGE_X_ERR ** 2, libcleoScan.X_Y_CORRELATION * libcleoScan.WORLD_AVERAGE_X_ERR * libcleoScan.WORLD_AVERAGE_Y_ERR),
           (libcleoScan.X_Y_CORRELATION * libcleoScan.WORLD_AVERAGE_X_ERR * libcleoScan.WORLD_AVERAGE_Y_ERR, libcleoScan.WORLD_AVERAGE_Y_ERR ** 2))

    vals = np.random.multivariate_normal(mean, cov, size=n)
    return vals[:, 0], vals[:, 1]


def _plot_xy(x_vals, y_vals, path):
    plt.hist2d(x_vals, y_vals, bins=50)
    plt.xlabel("x")
    plt.xlabel("y")
    plt.colorbar()
    plt.savefig(path)
    plt.clf()


def _plot_fit(actual_params, fit_params, ratios, errors, bin_centres, bin_widths, path):
    plt.clf()

    a, b, c = libcleoScan.expectedParams(actual_params)
    a_fit, b_fit, c_fit = libcleoScan.expectedParams(fit_params)

    plt.errorbar(bin_centres, ratios, yerr=errors, xerr=bin_widths/2, fmt="b+")
    plt.plot(bin_centres, [a_fit + b_fit*x + c_fit*x*x for x in bin_centres], "r", label="Fit")
    plt.plot(bin_centres, [a + b*x + c*x*x for x in bin_centres], "k--", label="Ideal")
    plt.legend()
    plt.xlabel("time")
    plt.ylabel(r"$\frac{WS}{RS}$ ratio")
    plt.savefig(path)


def _fit(i,
        k,
        actual_r,
        actual_re_z,
        x,
        y,
        width,
        bins,
        max_time,
        n_rs,
        out_q):

    # No efficiency
    efficiency = lambda x: 1

    decay_params = [x, y, actual_r, 0.7, actual_re_z, width]

    # Create a fitter
    Fitter = libcleoScan.ConstrainedFitter(bins, decay_params, [1]*6)

    for _ in range(k):
         seed = np.random.randint(0, np.iinfo(np.uint32).max)
         rs_times, ws_times = libcleoScan.simulate(n_rs, decay_params, max_time, seed)

         # Bin
         Fitter.addRSPoints(rs_times, np.ones(len(rs_times)))
         Fitter.addWSPoints(ws_times, np.ones(len(ws_times)))

    # Perform the fit
    Fitter.fixParameters(["z_im", "width"])
    result = Fitter.fit(efficiency)

    # _plot_fit(decay_params,
    #           result.fitParams,
    #           np.array(Fitter.ratios()),
    #           np.array(Fitter.errors()),
    #           np.array(Fitter.getBinCentres()),
    #           np.array(Fitter.getBinWidths()),
    #           f"{i}.png")
    # print(i)
    # sys.exit(0)
    out_dict = {i: (result.fitParams[2], result.fitParamErrors[2], result.fitParams[4], result.fitParamErrors[4])}

    out_q.put(out_dict)


def _plot_pulls(re_z_vals, r_vals, path):
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    kw = {"histtype": "step", "bins": np.linspace(-3, 3)}

    ax[0].hist(re_z_vals, **kw)
    ax[0].set_title(fr"{np.mean(re_z_vals):4.3}$\pm${np.std(re_z_vals):4.3}")
    ax[0].set_xlabel(r"$Re(Z)$")

    ax[0].text(0.5, .5, f"{len(re_z_vals)}")

    ax[1].hist(r_vals, **kw)
    ax[1].set_title(fr"{np.mean(r_vals):4.3}$\pm${np.std(r_vals):4.3}")
    ax[1].set_xlabel(r"$r$")

    ax[1].text(0.5, 0.5, f"{len(r_vals)}")

    plt.savefig(path)
    plt.clf()


def _multiproc_pull_study(num_experiments,
                          x_vals,
                          y_vals,
                          k,
                          actual_r,
                          actual_re_z,
                          width,
                          bins,
                          max_time,
                          n_rs):

    i = 0
    chunk_size = 8
    resultdict = {}
    while i < num_experiments:
        out_q = Queue()
        procs = []

        num_in_q = 0
        for _ in range(chunk_size):
            p = Process(target=_fit, args=(i, k, actual_r, actual_re_z, x_vals[i], y_vals[i], width, bins, max_time, n_rs, out_q))
            procs.append(p)
            p.start()
            num_in_q += 1

            i += 1
            if i == num_experiments:
                break

        for _ in range(num_in_q):
            resultdict.update(out_q.get())

        for p in procs:
            p.join()
        print(i)

    r_vals = np.array([x[0] for x in resultdict.values()])
    r_errs = np.array([x[1] for x in resultdict.values()])
    re_z_vals = np.array([x[2] for x in resultdict.values()])
    re_z_errs = np.array([x[3] for x in resultdict.values()])

    re_z_vals = (re_z_vals - actual_re_z) / re_z_errs
    r_vals = (r_vals - actual_r) / r_errs

    return r_vals, r_errs, re_z_vals, re_z_errs


def main():
    # Define our parameters
    width = 2.5
    max_time = 5.0 / width
    n_rs = 10_000_000  # Number of RS evts to generate each time
    k = 50  # Number of times to repeat generation
    num_experiments = 1000
    actual_re_z, actual_r = -0.1, 0.055

    # Define bins
    bins = np.logspace(np.log10(0.01), np.log10(max_time), 5)
    bins[0] = 0.0
    bins = np.linspace(0, max_time, 15)

    # Find vals of x, y to use
    x_vals, y_vals = _gen_x_y(num_experiments)

    r_vals, r_errs, re_z_vals, re_z_errs = _multiproc_pull_study(num_experiments,
                                                                 x_vals,
                                                                 y_vals,
                                                                 k,
                                                                 actual_r,
                                                                 actual_re_z,
                                                                 width,
                                                                 bins,
                                                                 max_time,
                                                                 n_rs)
    print(r_vals, r_errs, re_z_vals, re_z_errs)

    print(f"{np.mean(re_z_vals)=}\t{np.std(re_z_vals)=}")
    print(f"{np.mean(r_vals)=}\t{np.std(r_vals)=}")

    _plot_pulls(re_z_vals, r_vals, "pull.png")


if __name__ == "__main__":
    main()

