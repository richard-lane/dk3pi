"""
Plot a projection of phsp before/after reweighting

For my STFC summer school poster

Uses the same data for train+test since it's just illustrative

"""

import sys, os
import matplotlib.pyplot as plt
import numpy as np
import script_util

sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
import reweight_utils
import reweighting


def main():
    # Read flat phsp
    file_name = "2018MCflat.root"
    tree_name = "TupleDstToD0pi_D0ToKpipipi/DecayTree"
    mc = reweight_utils.read_invariant_masses(
        file_name,
        tree_name,
        ("K_PX", "K_PY", "K_PZ", "K_PE"),
        ("pi1_PX", "pi1_PY", "pi1_PZ", "pi1_PE"),
        ("pi2_PX", "pi2_PY", "pi2_PZ", "pi2_PE"),
        ("pi3_PX", "pi3_PY", "pi3_PZ", "pi3_PE"),
    )

    # Perfom cuts on ipchi2 and stuff
    dst_bkgcat = reweight_utils.read_branch(file_name, tree_name, "Dstar_BKGCAT")
    d_bkgcat = reweight_utils.read_branch(file_name, tree_name, "D_BKGCAT")
    to_remove = []
    for i, (dst, d) in enumerate(zip(dst_bkgcat, d_bkgcat)):
        if dst or d:
            to_remove.append(i)
    mc = np.delete(mc, to_remove, axis=0)

    # Generate some flat phsp with lots of events
    k, pi1, pi2, pi3 = script_util.flat_phsp_points(250000)
    for i in range(len(k[0])):
        pi1.T[i], pi2.T[i] = reweight_utils.momentum_order(k.T[i], pi1.T[i], pi2.T[i])
    flat = reweight_utils.invariant_mass_parametrisation(k, pi1, pi2, pi3)

    # Train BDT
    print("training bdt")
    bdt = reweighting.init(flat, mc, n_estimators=150)

    print("reweighting")
    weights = reweighting.predicted_weights(bdt, mc)

    # Plot a projection of one of the invariant masses
    kw = {"bins": 150, "alpha": 0.3, "density": True}
    plt.hist(mc[:, 3], **kw, label="LHCb MC", color="blue")
    plt.hist(flat[:, 3], **kw, label="Phsp", color="darkorange")
    plt.yticks([])
    plt.xlabel(r"$M(K\pi_1\pi_2) /MeV$")
    plt.legend()
    plt.show()

    # Plot the reweighted projection
    plt.hist(mc[:, 3], **kw, weights=weights, label="Reweighted LHCb MC", color="blue")
    plt.hist(flat[:, 3], **kw, label="Phsp", color="darkorange")
    plt.yticks([])
    plt.xlabel(r"$M(K\pi_1\pi_2) /MeV$")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
