"""
Create a scatter plot of random numbers

"""
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import uproot

plt.style.use("fivethirtyeight")
sys.path.append(
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "AmpGenTools",
        "python",
    )
)
sys.path.append(
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "efficiency",
        "scripts",
    )
)
import amplitude_models
import script_util


def _polar_plot(magnitudes, phases):
    """
    Plot a polar scatter plot of mag/phase, using some bin limits

    deprecated but keeping it in just in case

    """
    bin_limits = (-180.0, -39.0, 0.0, 43.0, 180.0)
    binned_phases = [[], [], [], []]
    binned_mags = [[], [], [], []]

    # Bin them into distinct bins
    for mag, phase in zip(magnitudes, phases):
        for i in range(4):
            if bin_limits[i] < phase <= bin_limits[i + 1]:
                binned_phases[i].append(phase)
                binned_mags[i].append(mag)
                break

    # Plot them on a polar plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="polar")
    for mags, phases in zip(binned_mags, binned_phases):
        ax.plot(np.deg2rad(phases), mags, ".", markersize=1)

    for lim in bin_limits[:-1]:
        ax.plot((0, np.deg2rad(lim)), (0, 1), "k")

    ax.set_ylim((0, 1))
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.grid(False)
    ax.spines["polar"].set_visible(False)
    plt.show()


def _bin_phases(phases):
    bin_limits = (-180.0, -39.0, 0.0, 43.0, 180.0)

    binned_phases = [[], [], [], []]
    for phase in phases:
        for i in range(len(bin_limits) - 1):
            if bin_limits[i] < phase <= bin_limits[i + 1]:
                binned_phases[i].append(phase)
                break

    return binned_phases


def _hist(phases, title):
    """
    Phases in deg

    """
    binned_phases = _bin_phases(phases)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for phases in binned_phases:
        ax.hist(np.deg2rad(phases), bins=np.linspace(-np.pi, np.pi, 100))

    ax.set_yticklabels([])
    ax.grid(False)
    ax.set_xlabel(r"$\delta$")

    x_tick = np.array((-1, -0.5, 0, 0.5, 1))
    ax.set_xticks(x_tick * np.pi)
    ax.set_xticklabels([fr"${x}\pi$" for x in x_tick])

    plt.title(title)
    plt.savefig(title.replace(" ", "_") + ".png", bbox_inches="tight")


def _read(file_name):
    """
    """
    tree = uproot.open(file_name)["DalitzEventList"]
    return np.row_stack(
        (
            tree.array("_1_K~_Px"),
            tree.array("_1_K~_Py"),
            tree.array("_1_K~_Pz"),
            tree.array("_1_K~_E"),
            tree.array("_2_pi#_Px"),
            tree.array("_2_pi#_Py"),
            tree.array("_2_pi#_Pz"),
            tree.array("_2_pi#_E"),
            tree.array("_3_pi#_Px"),
            tree.array("_3_pi#_Py"),
            tree.array("_3_pi#_Pz"),
            tree.array("_3_pi#_E"),
            tree.array("_4_pi~_Px"),
            tree.array("_4_pi~_Py"),
            tree.array("_4_pi~_Pz"),
            tree.array("_4_pi~_E"),
        )
    )


def main():
    # Read in RS, WS, flat data (use ampgen)
    rs = _read("Favoured.root")
    phases = [np.rad2deg(amplitude_models.phase(evt, 1)) for evt in rs.T]
    _hist(phases, "AmpGen RS phases")

    ws = _read("Mixed.root")
    phases = [np.rad2deg(amplitude_models.phase(evt, 1)) for evt in ws.T]
    _hist(phases, "AmpGen WS phases")

    flat = np.row_stack(script_util.flat_phsp_points(1000000)) / 1000.0
    phases = [np.rad2deg(amplitude_models.phase(evt, 1)) for evt in flat.T]
    _hist(phases, "Phsp Evts phases")

    # _polar_plot(magnitudes, phases)


if __name__ == "__main__":
    main()
