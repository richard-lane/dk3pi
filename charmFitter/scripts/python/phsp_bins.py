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
import reweight_utils


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


def _binned_invmasses(events, phases, *particles):
    """
    Find the invariant masses of a number of particles

    particles should be any of 0,1,2,3

    """
    bin_limits = (-180.0, -39.0, 0.0, 43.0, 180.0)
    x = [[] for _ in range(4)]

    for evt, phase in zip(events.T, phases):
        px = 0.0
        py = 0.0
        pz = 0.0
        e = 0.0
        for particle in particles:
            px += evt[0 + particle * 4]
            py += evt[1 + particle * 4]
            pz += evt[2 + particle * 4]
            e += evt[3 + particle * 4]

        s = reweight_utils.invariant_mass(px, py, pz, e)

        # Append to the right list based on phase
        for i in range(len(bin_limits) - 1):
            if bin_limits[i] < phase <= bin_limits[i + 1]:
                x[i].append(s)

    return x


def _dalitz(
    rs, ws, flat, rs_phases, ws_phases, flat_phases, title, x_indices, y_indices
):
    """
    Plot a scatter plot of points, colour coded by phsp bin

    finds invariant masses s_ij based on ij... in x_indices, y_indices

    """
    # Plot
    fig, ax = plt.subplots(1, 3)

    for a in ax:
        a.set_xlabel(r"$M(K\pi_1)$")
        a.grid(False)
        a.set_xticklabels([])
        a.set_yticklabels([])

    ax[0].set_ylabel(r"$M(K\pi_2\pi_3)$")

    for i, (data, phase) in enumerate(
        zip((rs, ws, flat), (rs_phases, ws_phases, flat_phases))
    ):
        binned_data_x = _binned_invmasses(data, phase, *x_indices)
        binned_data_y = _binned_invmasses(data, phase, *y_indices)
        for (dx, dy) in zip(binned_data_x, binned_data_y):
            ax[i].plot(dx, dy, ".")

    fig.suptitle(title)
    plt.show()


def _read(file_name):
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


def _all_two_body_projections(rs, ws, flat, rs_phases, ws_phases, flat_phases):
    # lol
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    # Plot only interesting ones
                    if a != b and c != d:
                        if set((a, b)) != set((c, d)):
                            _dalitz(
                                rs,
                                ws,
                                flat,
                                rs_phases,
                                ws_phases,
                                flat_phases,
                                r"$s_{"
                                + str(a)
                                + str(b)
                                + r"} vs s_{"
                                + str(c)
                                + str(d)
                                + r"}$",
                                (a, b),
                                (c, d),
                            )


def _simple_binning_plot(data):
    x = [[] for _ in range(4)]
    y = [[] for _ in range(4)]

    lims = 1.0, 1.4

    # Find inv masses
    for evt in data.T:
        s01 = reweight_utils.invariant_mass(
            evt[0] + evt[4], evt[1] + evt[5], evt[2] + evt[6], evt[3] + evt[7]
        )
        s023 = reweight_utils.invariant_mass(
            evt[0] + evt[8] + evt[12],
            evt[1] + evt[9] + evt[13],
            evt[2] + evt[10] + evt[14],
            evt[3] + evt[11] + evt[15],
        )
        bin = 0
        if s01 < lims[0] and s023 > lims[1]:
            bin=1
        if s01 > lims[0] and s023 < lims[1]:
            bin=2
        if s01 > lims[0] and s023 > lims[1]:
            bin=3

        x[bin].append(s01)
        y[bin].append(s023)

    # Plot
    for i in range(4):
        plt.plot(x[i], y[i], ".")

    plt.xlabel(r"$M(K\pi_1) /GeV$")
    plt.ylabel(r"$M(K\pi_2\pi_3) /GeV$")
    plt.savefig("phsp.png", bbox_inches="tight")


def main():
    # Read in RS, WS, flat data (use ampgen)
    #rs = _read("shortFavoured.root")
    #rs_phases = [np.rad2deg(amplitude_models.phase(evt, 1)) for evt in rs.T]

    #ws = _read("shortMixed.root")
    #ws_phases = [np.rad2deg(amplitude_models.phase(evt, 1)) for evt in ws.T]

    flat = np.row_stack(script_util.flat_phsp_points(3000)) / 1000.0
    print(flat)
    #flat_phases = [np.rad2deg(amplitude_models.phase(evt, 1)) for evt in flat.T]

    _simple_binning_plot(flat)

    # Two body dalitz plots
    # _all_two_body_projections(rs, ws, flat, rs_phases, ws_phases, flat_phases)


if __name__ == "__main__":
    main()
