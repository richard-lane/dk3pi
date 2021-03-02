"""
Create a scatter plot of random numbers

"""
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("fivethirtyeight")


def main():
    # Create a load of random phases and magnitudes
    N = 50000
    magnitudes = np.random.triangular(0, 1.0, 1.0, N)
    phases = 360.0 * np.random.random(N) - 180.0

    # Bin them into distinct bins
    bin_limits = (-180.0, -39.0, 0.0, 43.0, 180.0)

    binned_phases = [[], [], [], []]
    binned_mags = [[], [], [], []]

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
    #ax.spines["polar"].set_visible(False)
    plt.show()


if __name__ == "__main__":
    main()
