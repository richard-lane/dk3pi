"""
Helper scripts for common plots

"""
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs


def plot_efficiencies(domain, *efficiencies, labels=None, path=None):
    """
    Plot a collection of efficiency projections

    :param domain:       5-member iterable of domains over which the efficincies are defined. Likely bin centres
    :param efficiencies: Iterable of 5 efficiency projections and their errors.
                         Each should be a tuple of 5 array-likes of efficiencies (Kpi1, pi1pi2, pi2pi3, kpi1pi2, pi1pi2pi3) and their errors,
                         i.e. should have shape (2, 5)
    :param labels:       A label for each efficiency, if specified.
    :param path:         Path to save to if specified. If not, will display the plot instead.

    """
    dim = 5
    assert len(domain) == dim, f"{dim} members in domain"
    for e in efficiencies:
        assert len(e[0]) == dim, f"{dim} iterables of efficiencies"
        assert len(e[1]) == dim, f"{dim} iterables of errors"

    plt.figure(figsize=[12.8, 9.6])
    gs.GridSpec(2, 3)

    # First plot as title
    plt.subplot2grid((2, 3), (0, 0))
    plt.text(0.2, 0.5, "Efficiencies", fontsize=16, bbox={"color": "white"})
    plt.axis("off")

    for i, ax, ax_label in zip(
        range(dim),
        ([0, 1], [0, 2], [1, 0], [1, 1], [1, 2]),
        (
            r"M($K\pi_1$)",
            r"M($\pi_1\pi_2$)",
            r"M($\pi_2\pi_3$)",
            r"M($K\pi_1\pi_2$)",
            r"M($\pi_1\pi_2\pi_3$)",
            r" ",
        ),
    ):
        plt.subplot2grid((2, 3), ax)
        for j, e in enumerate(efficiencies):
            label = None if labels is None else labels[j]
            plt.plot(domain[i], e[0][i], label=label)
            plt.fill_between(domain[i], e[0][i] + e[1][i], e[0][i] - e[1][i], alpha=0.5)
        plt.title(ax_label)
        plt.ylim((0, 2))

        # Only need one legend
        if not i:
            plt.legend()

        # Bottom two need axis labels
        if i in {3, 4}:
            plt.xlabel("MeV")

    if path:
        plt.savefig(path, dpi=600)
    else:
        plt.show()


def classification_plots(
    mc,
    model,
    mc_prob,
    model_prob,
    title,
    label,
    chisq,
    p,
    prob_bins,
    mc_weights=None,
    path=None,
):
    """
    Make plots of the phsp projections of our flat and MC data, and also plot the
    classification probabilities

    TODO document this haha lol

    """
    plt.figure()
    gs.GridSpec(3, 4)
    plt.suptitle(title)

    # Use the first plot as a title
    plt.subplot2grid((3, 4), (0, 0))
    plt.text(0.2, 0.5, "Phase Space\nProjections", fontsize=16, bbox={"color": "white"})
    plt.axis("off")

    # Large plot to show classification probabilities
    kw = {"alpha": 0.3}
    plt.subplot2grid((3, 4), (0, 2), colspan=2, rowspan=3)
    plt.hist(mc_prob, **kw, label="MC", bins=prob_bins, weights=mc_weights)
    plt.hist(model_prob, **kw, label=label, bins=prob_bins)
    plt.legend()
    plt.yticks([])
    plt.title("Classification probability")
    plt.text(0.1, 0.8, f"$\chi^2$={chisq}\np={p}", transform=plt.gca().transAxes)

    # Small plots to show phsp projections
    for i, ax, ax_label in zip(
        range(5),
        ([0, 1], [1, 0], [1, 1], [2, 0], [2, 1]),
        (
            r"M($K\pi_1$)",
            r"M($\pi_1\pi_2$)",
            r"M($\pi_2\pi_3$)",
            r"M($K\pi_1\pi_2$)",
            r"M($\pi_1\pi_2\pi_3$)",
            r" ",
        ),
    ):
        plt.subplot2grid((3, 4), ax)
        plt.hist(mc[:, i], **kw, label="MC", weights=mc_weights, bins=75)
        plt.hist(model[:, i], **kw, label=label, bins=75)
        plt.title(ax_label)

        # Only need one legend
        if not i:
            plt.legend()

        # Bottom two need axis labels
        if i in {3, 4}:
            plt.xlabel("MeV")

    if path:
        plt.savefig(path, dpi=600)
    else:
        plt.show()
