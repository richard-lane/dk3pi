"""
Pull study for fitter

"""
# Set up include dir
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from multiprocessing import Queue, Process
from scipy.stats import norm

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), "build/charmFitter/python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), "dk3pi-efficiency-model/"))

import libcleoScan
from model.util import definitions


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
    Fitter = libcleoScan.ConstrainedFitter(bins, decay_params, [1, 1, 1, 1, 1, 1])

    for _ in range(k):
        # Reset numpy's seed - do this explicitly so different processes will
        # generate different random numbers
        np.random.seed()

        seed = np.random.randint(0, np.iinfo(np.uint32).max)
        rs_times, ws_times = libcleoScan.simulate(n_rs, decay_params, max_time, seed)

        # Bin
        Fitter.addRSPoints(rs_times, np.ones(len(rs_times)))
        Fitter.addWSPoints(ws_times, np.ones(len(ws_times)))

    # Perform the fit with fixed Re(Z)
    Fitter.fixParameters(["z_re", "z_im", "width"])
    result_fixed_rez = Fitter.fit(efficiency)

    # Perform the fit again with fixed rD
    Fitter.freeParameter("z_re")
    Fitter.setParameter("z_re", actual_re_z)
    Fitter.fixParameter("r")
    result_fixed_r = Fitter.fit(efficiency)

    # Perform again with free Re(Z) and rD
    Fitter.freeParameter("r")
    Fitter.setParameter("r", actual_r)
    Fitter.setParameter("z_re", actual_re_z)
    Fitter.setParameter("x", x)
    Fitter.setParameter("y", y)
    result = Fitter.fit(efficiency)

    rez_chi2_diff = np.sqrt(abs(result.fitStatistic - result_fixed_rez.fitStatistic))
    r_chi2_diff = np.sqrt(abs(result.fitStatistic - result_fixed_r.fitStatistic))

    # _plot_fit(decay_params,
    #           result.fitParams,
    #           np.array(Fitter.ratios()),
    #           np.array(Fitter.errors()),
    #           np.array(Fitter.getBinCentres()),
    #           np.array(Fitter.getBinWidths()),
    #           f"{i}.png")
    # print(i)
    # sys.exit(0)
    out_dict = {i: (result.fitParams[2], result.fitParamErrors[2], result.fitParams[4], result.fitParamErrors[4], rez_chi2_diff, r_chi2_diff)}

    out_q.put(out_dict)


def _plot_pulls(re_z_vals, r_vals, path):
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    kw = {"histtype": "step", "bins": np.linspace(-3, 3), "density": True}
    centres = (kw["bins"][1:] + kw["bins"][:-1]) / 2.0

    ax[0].hist(re_z_vals, **kw)
    # Overlay ideal + actual bell curves
    mean, std = np.mean(re_z_vals), np.std(re_z_vals)
    ax[0].plot(centres, [norm.pdf(x) for x in centres], "k--", alpha=0.6)
    ax[0].plot(centres, [norm.pdf(x, mean, std) for x in centres], "r")

    ax[0].set_title(fr"Re(Z): {mean:4.3}$\pm${std:4.3}")
    ax[0].set_xlabel(r"$\frac{x-\mu}{\sigma}$")

    ax[1].hist(r_vals, **kw)

    mean, std = np.mean(r_vals), np.std(r_vals)
    ax[1].plot(centres, [norm.pdf(x) for x in centres], "k--", alpha=0.6)
    ax[1].plot(centres, [norm.pdf(x, mean, std) for x in centres], "r")

    ax[1].set_title(fr"$r_D: ${mean:4.3}$\pm${std:4.3}")
    ax[1].set_xlabel(r"$\frac{x-\mu}{\sigma}$")

    ax[1].text(0.5, 0.5, f"num experiments: {len(r_vals)}")

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
    rez_chi2_diffs = np.array([x[4] for x in resultdict.values()])
    r_chi2_diffs = np.array([x[5] for x in resultdict.values()])

    re_z_vals = (re_z_vals - actual_re_z) / re_z_errs
    r_vals = (r_vals - actual_r) / r_errs

    return r_vals, r_errs, re_z_vals, re_z_errs, rez_chi2_diffs, r_chi2_diffs


def _plot_coverage(chi2_diffs, label, title, path):
    contents, bins, _ = plt.hist(chi2_diffs, cumulative=True, bins=100, histtype="step", label=label)
    centres = (bins[1:] + bins[:-1] ) / 2.0

    # Normal CDF is 0.5 at 0
    expected = [2 * len(chi2_diffs) * (norm.cdf(x) - 0.5) for x in centres]
    plt.plot(centres, expected, "k+", label="Gaussian")

    plt.title(title)
    plt.xlabel(r"$\Delta \chi^2$")
    plt.ylabel("num experiments")

    plt.legend()

    print(f"Saving {path}")
    plt.savefig(path)
    plt.clf()


def main():
    # Define our parameters
    width = 1.0 / definitions.D_LIFETIME
    max_time = definitions.D_LIFETIME * definitions.TIME_BINS[-1]
    # We expect ~60M RS evts in each phsp bin...
    n_rs = 1_000_000  # Number of RS evts to generate each time
    k = 6  # Number of times to repeat generation
    num_experiments = 500
    actual_re_z, actual_r = -0.1, 0.055

    # Define bins
    bins = [definitions.D_LIFETIME * x for x in definitions.TIME_BINS[1:]]

    # Find vals of x, y to use
    x_vals, y_vals = _gen_x_y(num_experiments)

    r_vals, r_errs, re_z_vals, re_z_errs, rez_chi2_diffs, r_chi2_diffs = _multiproc_pull_study(num_experiments,
                                                                 x_vals,
                                                                 y_vals,
                                                                 k,
                                                                 actual_r,
                                                                 actual_re_z,
                                                                 width,
                                                                 bins,
                                                                 max_time,
                                                                 n_rs)
    print(f"{np.mean(re_z_vals)=}\t{np.std(re_z_vals)=}")
    print(f"{np.mean(r_vals)=}\t{np.std(r_vals)=}")

    _plot_pulls(re_z_vals, r_vals, "pull.png")
    _plot_coverage(rez_chi2_diffs, r"$\Delta \chi^2$", r"$\Delta \chi^2$ between fits with fixed/free Re(Z)", "rez_coverage.png")
    _plot_coverage(r_chi2_diffs, r"$\Delta \chi^2$", r"$\Delta \chi^2$ between fits with fixed/free $r_D$", "r_coverage.png")


if __name__ == "__main__":
    main()

