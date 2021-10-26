"""
Simulate some RS and WS decays, create a fitter, do a fit, plot it

"""
# Set up include dir
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from scipy.integrate import quad

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "build/charmFitter/python"))

import libcleoScan


def main():
    # Define our parameters
    width = 2.5
    decay_params = [0.03, 0.06, 0.055, 0.5, 0.5, width]  # Set these to the world avg x, y
    max_time = 3.0
    n_rs = 10_000_000  # Number of RS evts to generate each time
    k = 50  # Number of times to repeat generation

    # Find ideal a, b, c
    a, b, c = libcleoScan.expectedParams(decay_params)

    # Define bins, find rates at bin centres
    bins = np.linspace(0, max_time, 25)

    # Create a fitter
    Fitter = libcleoScan.CLEOCombinationFitter(bins, decay_params, [1]*6, 0)

    # Simulate decays
    n_ws = 0
    for i in range(k):
         print("\tgenerating")
         seed = np.random.randint(0, np.iinfo(np.uint32).max)
         rs_times, ws_times = libcleoScan.simulate(n_rs, decay_params, max_time, seed)

         # Bin
         print("\tadding to fitter")
         Fitter.addRSPoints(rs_times, np.ones(len(rs_times)))
         Fitter.addWSPoints(ws_times, np.ones(len(ws_times)))

         print(f"done {i}")

    # Perform the fit
    print(f"{Fitter.getRSBinContent()=}")
    print(f"{Fitter.getWSBinContent()=}")
    print("fitting")
    Fitter.fixParameters(["z_im", "width", "r"])
    result = Fitter.fit(lambda x: 1)
    a_fit, b_fit, c_fit = libcleoScan.expectedParams(result.fitParams)

    # Get the fit results
    for i in range(6):
        print(f"{decay_params[i]}\t{result.fitParams[i]}+-{result.fitParamErrors[i]}\n")
    print(f"{a}\t{a_fit}")
    print(f"{b}\t{b_fit}")
    print(f"{c}\t{c_fit}")

    # Plot the resulting graph
    centres = np.array(Fitter.getBinCentres())
    widths = np.array(Fitter.getBinWidths())
    ratios = np.array(Fitter.ratios())
    errors = np.array(Fitter.errors())

    plt.errorbar(centres, ratios, yerr=errors, xerr=widths/2, fmt="b+")
    plt.plot(centres, [a_fit + b_fit*x + c_fit*x*x for x in centres], "r", label="Fit")
    plt.plot(centres, [a + b*x + c*x*x for x in centres], "k--", label="Ideal")
    plt.legend()
    plt.xlabel("time")
    plt.ylabel(r"$\frac{WS}{RS}$ ratio")
    plt.savefig("fit_example.png")


if __name__ == "__main__":
    main()

