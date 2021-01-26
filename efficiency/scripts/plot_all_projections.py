"""
Plot all two- and three- body invariant mass projections of a K3pi system

"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import phasespace
from multiprocessing import Process

import matplotlib

matplotlib.use("Agg")

sys.path.append(os.path.dirname(__file__) + "/../bdt_reweighting/")
import reweight_utils


def three_body_histograms(
    path, mc_k, mc_pi1, mc_pi2, mc_pi3, model_k, model_pi1, model_pi2, model_pi3
):
    """
    Take arrays of particles' kinematic data

    Save plots of three-body invariant masses

    Will save plots to path + some stuff

    """
    num_particles = 4
    mc_particles = (mc_k, mc_pi1, mc_pi2, mc_pi3)
    model_particles = (model_k, model_pi1, model_pi2, model_pi3)
    particle_names = (r"$K$", r"$\pi_1$", r"$\pi_2$", r"$\pi_3$")
    hist_kwargs = {"bins": 200, "alpha": 0.3, "density": True}

    for i in range(num_particles):
        for j in range(i + 1, num_particles):
            for k in range(j + 1, num_particles):
                mc_invariant_masses = reweight_utils.invariant_masses(
                    mc_particles[i][0] + mc_particles[j][0] + mc_particles[k][0],
                    mc_particles[i][1] + mc_particles[j][1] + mc_particles[k][1],
                    mc_particles[i][2] + mc_particles[j][2] + mc_particles[k][2],
                    mc_particles[i][3] + mc_particles[j][3] + mc_particles[k][3],
                )

                model_invariant_masses = reweight_utils.invariant_masses(
                    model_particles[i][0]
                    + model_particles[j][0]
                    + model_particles[k][0],
                    model_particles[i][1]
                    + model_particles[j][1]
                    + model_particles[k][1],
                    model_particles[i][2]
                    + model_particles[j][2]
                    + model_particles[k][2],
                    model_particles[i][3]
                    + model_particles[j][3]
                    + model_particles[k][3],
                )

                plt.hist(mc_invariant_masses, **hist_kwargs, label="MC")
                plt.hist(model_invariant_masses, **hist_kwargs, label="Model")
                plt.legend()
                label = f"{particle_names[i]}{particle_names[j]}{particle_names[k]}"
                plt.title(f"Inv mass: {label}")
                plt.xlabel(f"M({label}) /MeV")

                plt.savefig(f"{path}{i}{j}{k}.png")
                plt.clf()


def two_body_histograms(
    path, mc_k, mc_pi1, mc_pi2, mc_pi3, model_k, model_pi1, model_pi2, model_pi3
):
    """
    Take arrays of particles' kinematic data

    Save plots of two-body invariant masses

    """
    num_particles = 4
    mc_particles = (mc_k, mc_pi1, mc_pi2, mc_pi3)
    model_particles = (model_k, model_pi1, model_pi2, model_pi3)
    particle_names = (r"$K$", r"$\pi_1$", r"$\pi_2$", r"$\pi_3$")
    hist_kwargs = {"bins": 200, "alpha": 0.3, "density": True}

    for i in range(num_particles):
        for j in range(i + 1, num_particles):
            mc_invariant_masses = reweight_utils.invariant_masses(
                mc_particles[i][0] + mc_particles[j][0],
                mc_particles[i][1] + mc_particles[j][1],
                mc_particles[i][2] + mc_particles[j][2],
                mc_particles[i][3] + mc_particles[j][3],
            )
            model_invariant_masses = reweight_utils.invariant_masses(
                model_particles[i][0] + model_particles[j][0],
                model_particles[i][1] + model_particles[j][1],
                model_particles[i][2] + model_particles[j][2],
                model_particles[i][3] + model_particles[j][3],
            )

            plt.hist(mc_invariant_masses, **hist_kwargs, label="MC")
            plt.hist(model_invariant_masses, **hist_kwargs, label="Model")
            plt.legend()
            label = f"{particle_names[i]}{particle_names[j]}"
            plt.title(f"Inv mass: {label}")
            plt.xlabel(f"M({label}) /MeV")

            plt.savefig(f"{path}{i}{j}.png")
            plt.clf()


def mc_particles(file_name):
    """
    Read particle data from a MC file

    """
    tree_name = "TupleDstToD0pi_D0ToKpipipi/DecayTree"
    branches = (
        ("K_PX", "K_PY", "K_PZ", "K_PE"),
        ("pi1_PX", "pi1_PY", "pi1_PZ", "pi1_PE"),
        ("pi2_PX", "pi2_PY", "pi2_PZ", "pi2_PE"),
        ("pi3_PX", "pi3_PY", "pi3_PZ", "pi3_PE"),
    )
    return reweight_utils.read_kinematic_data(file_name, tree_name, *branches)


def ampgen_particles(file_name, label):
    """
    Read particle data from an ampgen file

    """
    assert label in ("rs", "ws")

    tree_name = "DalitzEventList"
    branches = (
        (
            ("_1_K#_Px", "_1_K#_Py", "_1_K#_Pz", "_1_K#_E"),
            ("_2_pi~_Px", "_2_pi~_Py", "_2_pi~_Pz", "_2_pi~_E"),
            ("_3_pi~_Px", "_3_pi~_Py", "_3_pi~_Pz", "_3_pi~_E"),
            ("_4_pi#_Px", "_4_pi#_Py", "_4_pi#_Pz", "_4_pi#_E"),
        )
        if label == "rs"
        else (
            ("_1_K~_Px", "_1_K~_Py", "_1_K~_Pz", "_1_K~_E"),
            ("_2_pi#_Px", "_2_pi#_Py", "_2_pi#_Pz", "_2_pi#_E"),
            ("_3_pi#_Px", "_3_pi#_Py", "_3_pi#_Pz", "_3_pi#_E"),
            ("_4_pi~_Px", "_4_pi~_Py", "_4_pi~_Pz", "_4_pi~_E"),
        )
    )
    # Convert to MeV by *1000
    return [
        1000 * x
        for x in reweight_utils.read_kinematic_data(file_name, tree_name, *branches)
    ]


def rs():
    mc = mc_particles("2018MC_RS.root")
    ampgen = ampgen_particles("ampgen_RS.root", "rs")
    two_body_histograms("rs/plot", *mc, *ampgen)
    three_body_histograms("rs/plot", *mc, *ampgen)


def ws():
    mc = mc_particles("2018MC_WS.root")
    ampgen = ampgen_particles("ampgen_WS.root", "ws")
    two_body_histograms("ws/plot", *mc, *ampgen)
    three_body_histograms("ws/plot", *mc, *ampgen)


def rs_dbar():
    """
    Compare RS LHCb MC with the AmpGen RS Dbar model

    """
    mc = mc_particles("2018MC_RS.root")
    ampgen = ampgen_particles("ampgen_Dbar_RS.root", "ws")
    two_body_histograms("rs_dbar/plot", *mc, *ampgen)
    three_body_histograms("rs_dbar/plot", *mc, *ampgen)


def flat():
    # Read phsp mc data
    mc = mc_particles("2018MCflat.root")

    # Generate flat data
    pi_mass = 139.570
    k_mass = 493.677
    d_mass = 1864.84

    # Find the kinematic information
    # Could possibly make this more efficient by doing it in a loop, maybe
    num_events = 500000
    _, particles = phasespace.nbody_decay(
        d_mass, (k_mass, pi_mass, pi_mass, pi_mass)
    ).generate(n_events=num_events)

    # Create labels for the particles to make it easier to code because im tired
    flat = (
        particles["p_0"].numpy().T,
        particles["p_1"].numpy().T,
        particles["p_2"].numpy().T,
        particles["p_3"].numpy().T,
    )

    two_body_histograms("flat/plot", *mc, *flat)
    three_body_histograms("flat/plot", *mc, *flat)


def main():
    flat_plots = Process(target=flat)
    rs_plots = Process(target=rs)
    ws_plots = Process(target=ws)
    rs_dbar_plots = Process(target=rs_dbar)

    flat_plots.start()
    rs_plots.start()
    ws_plots.start()
    rs_dbar_plots.start()

    flat_plots.join()
    rs_plots.join()
    ws_plots.join()
    rs_dbar_plots.join()


if __name__ == "__main__":
    main()
