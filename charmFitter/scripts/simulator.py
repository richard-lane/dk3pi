"""
Simulate RS/WS D->K3pi decay times

"""
# Set up include dir
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from scipy.integrate import quad

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "build/charmFitter/python"))

import libcleoScan


def normalised_rates(max_time, width, a, b, c):
    """
    Normalised RS/WS decay rates; i.e. PDFs

    """
    raw_rs_rate = lambda t: np.exp(-width * t)
    raw_ws_rate = lambda t : (a + b*t + c*t*t) * raw_rs_rate(t)

    rs_integral = quad(raw_rs_rate, 0, max_time)[0]
    ws_integral = quad(raw_ws_rate, 0, max_time)[0]

    return lambda t: raw_rs_rate(t) / rs_integral, lambda t: raw_ws_rate(t) / ws_integral


def _hist(normalised_rate, counts, bins, title, path):
    plt.clf()

    centres = (bins[1:] + bins[:-1]) / 2
    widths = bins[1:] - bins[:-1]

    n_bins = len(bins) - 1
    expected = widths * np.sum(counts) * np.array([normalised_rate(x) for x in centres])

    plt.hist(bins[:-1], bins, weights=counts, histtype="step", label="Simulated")
    plt.plot(centres, expected, "k--", linewidth=1, label="Expected")

    plt.title(title)
    plt.legend()
    plt.xlabel("time")
    plt.ylabel("counts")

    plt.savefig(path)


def _ratio(actual_ratios, ratio_errs, expected_ratios, bins, path):
    plt.clf()

    widths = bins[1:] - bins[:-1]
    centres = (bins[1:] + bins[:-1]) / 2

    plt.errorbar(centres, actual_ratios, xerr=widths, yerr=ratio_errs, fmt="b+", label="Simulated", alpha=0.5)
    plt.plot(centres, expected_ratios, "k--", label="Expected")

    plt.legend()
    plt.ylabel("WS/RS ratio")
    plt.xlabel("time")

    plt.savefig(path)


def _pdfs(rs_pdf, ws_pdf, bins, max_time, path):
    plt.clf()

    widths = bins[1:] - bins[:-1]
    centres = (bins[1:] + bins[:-1]) / 2
    n_bins = len(bins) - 1

    # Find pdf describing ratio
    ratio_pdf = lambda t: ws_pdf(t) / rs_pdf(t)
    integral = quad(ratio_pdf, 0, max_time)[0]
    normalised_ratio_pdf = lambda t: ratio_pdf(t) / integral

    # Plot RS
    naive_rs = widths * np.array([rs_pdf(x) for x in centres])
    exact_rs = np.array([quad(rs_pdf, bins[i], bins[i+1])[0] for i in range(n_bins)])

    # Plot WS
    naive_ws = widths * np.array([ws_pdf(x) for x in centres])
    exact_ws = np.array([quad(ws_pdf, bins[i], bins[i+1])[0] for i in range(n_bins)])

    # Plot ratio
    naive_ratio = widths * np.array([normalised_ratio_pdf(x) for x in centres])
    exact_ratio = np.array([quad(normalised_ratio_pdf, bins[i], bins[i+1])[0] for i in range(n_bins)])

    # Plot ratio of naive and exact
    naive = naive_ws / naive_rs
    exact = exact_ws / exact_rs

    fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    ax[0, 0].plot(centres, naive_rs, "r--", label="Approx Occupancy")
    ax[0, 0].plot(centres, exact_rs, "k.", label="Exact Occupancy", markersize=1.5)
    ax[0, 0].set_title("RS pdf")

    ax[0, 1].plot(centres, naive_ws, "r--", label="Approx Occupancy")
    ax[0, 1].plot(centres, exact_ws, "k.", label="Exact Occupancy", markersize=1.5)
    ax[0, 1].set_title("WS")

    ax[1, 0].plot(centres, naive_ratio, "r--", label="Approx Occupancy")
    ax[1, 0].plot(centres, exact_ratio, "k.", label="Exact Occupancy", markersize=1.5)
    ax[1, 0].set_title("Ratio (pdf)")

    ax[1, 1].plot(centres, naive, "r--", label="Approx Occupancy")
    ax[1, 1].plot(centres, exact, "k.", label="Exact Occupancy", markersize=1.5)
    ax[1, 1].set_title("Ratio of RS/WS pdfs")

    for a in ax.flatten():
        a.legend()

    plt.savefig(path)


def main():
    # Define our parameters
    width = 2.5
    decay_params = [0.03, 0.06, 0.055, 0.5, 0.5, width]
    max_time = 3.0
    n_rs = 10_000_000  # Number of RS evts to generate each time
    k = 50  # Number of times to repeat generation

    # Find ideal a, b, c
    a, b, c = libcleoScan.expectedParams(decay_params)

    # RS, WS rates
    rs_rate, ws_rate = normalised_rates(max_time, width, a, b, c)

    # Define bins, find rates at bin centres
    bins = np.linspace(0, max_time, 25)
    centres = (bins[:-1] + bins[1:]) / 2.0

    # Simulate decays
    n_bins = len(bins) - 1
    rs_binned, ws_binned = np.zeros(n_bins), np.zeros(n_bins)
    n_ws = 0
    for i in range(k):
         seed = np.random.randint(0, np.iinfo(np.uint32).max)
         rs_times, ws_times = libcleoScan.simulate(n_rs, decay_params, max_time, seed)

         # Bin
         rs_binned_tmp, _ = np.histogram(rs_times, bins)
         ws_binned_tmp, _ = np.histogram(ws_times, bins)

         rs_binned += rs_binned_tmp
         ws_binned += ws_binned_tmp
         print(f"done {i}")
         n_ws += len(ws_times)
    n_rs = k * n_rs

    print(f"{n_rs=}")
    print(f"{n_ws=}")

    # Plot histograms
    _hist(rs_rate, rs_binned, bins, "RS", "rs.png")
    _hist(ws_rate, ws_binned, bins, "WS", "ws.png")
    
    # num expected approx = width * rate at centre
    binwidths = bins[1:] - bins[:-1]
    naive_expected_rs = n_rs * binwidths * np.array([rs_rate(t) for t in centres])
    naive_expected_ws = n_ws * binwidths * np.array([ws_rate(t) for t in centres])

    # Num expected = integral
    expected_rs = n_rs * np.array([quad(rs_rate, bins[i], bins[i+1])[0] for i in range(n_bins)])
    expected_ws = n_ws * np.array([quad(ws_rate, bins[i], bins[i+1])[0] for i in range(n_bins)])

    # Ratios
    naive_expected_ratio = naive_expected_ws / naive_expected_rs
    expected_ratio = expected_ws / expected_rs

    actual_ratios = ws_binned / rs_binned
    ratio_errs = actual_ratios * np.sqrt(1/ws_binned+ 1/rs_binned)

    _ratio(actual_ratios, ratio_errs, naive_expected_ratio, bins, "naive.png")
    _ratio(actual_ratios, ratio_errs, expected_ratio, bins, "actual.png")

    # Plot each with idealised curve on each
    _pdfs(rs_rate, ws_rate, bins, max_time, "pdfs.png")


if __name__ == "__main__":
    main()

