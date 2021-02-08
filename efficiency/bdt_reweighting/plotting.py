"""
Helper scripts for common plots

"""
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs


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
        plt.yticks([])
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
