"""
Do a check for CLEO scan using the python binding

"""
# Set up include dir
import numpy as np
import sys, os
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "build/charmFitter/python"))

import libcleoScan
import common


def main():
    # Simulate decays
    # Define our parameters
    width = 2.5
    decay_params = [0.03, 0.06, 0.055, 0.5, 0.5, width]
    max_time = 3.0
    n_rs = 1_000_000  # Number of RS evts to generate each time
    k = 5  # Number of times to repeat generation
    bins = np.linspace(0, max_time, 15)
    
    # Simulate decays
    seed = np.random.randint(0, np.iinfo(np.uint32).max)
    rs_t, ws_t= libcleoScan.simulate(n_rs, decay_params, max_time, seed)

    # Set up params
    n_re_z, n_im_z = 25, 25
    reZVals = np.linspace(-0.9, 0.9, n_re_z)
    imZVals = np.linspace(-0.9, 0.9, n_im_z)
    rsWeights    = np.ones(len(rs_t))
    wsWeights    = np.ones(len(ws_t))
    initialErrs  = [1, 1, 1, 1, 1, 1]

    # Create a fitter
    Fitter = libcleoScan.ConstrainedFitter(bins, decay_params, initialErrs)
    Fitter.addRSPoints(rs_t, np.ones_like(rs_t))
    Fitter.addWSPoints(ws_t, np.ones_like(ws_t))
    Fitter.fixParameters(("width", "z_im", "z_re"))

    l = []
    for r in reZVals:
        for i in imZVals:
            Fitter.setParameter("z_re", r)
            Fitter.setParameter("z_im", i)
            l.append(Fitter.fit(lambda x: 1).fitStatistic)
            common.reset_params(Fitter, decay_params)
    l = np.array(l)

    # Subtract off minimum
    l = 2 * np.log(np.array(l))  # Should be -2?
    l -= np.min(l)
    l[l > 3.0] = 0.0

    # Convert to 2d array
    l = l.reshape((n_re_z, n_im_z), order="C")

    # Plot something
    print(l)
    plt.contourf(reZVals, imZVals, l, 100)
    plt.savefig("py_scan.png")


if __name__ == "__main__":
    main()

