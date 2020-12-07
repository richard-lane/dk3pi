"""
Helpers for visualising our multidimensional phase space distributions

"""


class Slices:
    def __init__(self, num_slices, num_bins, bin_limits, plot_index, slice_index):
        """
        Init an instance for plotting histogram projections of the plot_index'th parameter,
        taking slices along the slice_index'th parameter

        """
        pass

    def add_point(point, weight=1):
        pass

    def add_points(points, weights=None):
        pass

    # Array of histograms representing each slice through phase space
    _slice_array = None

    # Bin limits for finding which histogram a point belongs in
    _slice_limits = None

    # Indices of variables that we're using to slice across, and plot projections in
    _plot_index = None
    _slice_index = None

    # Track the total number of points in this slice
    # Float since the points are weighted
    _num_points = 0.0


def plot_slices(path, slices, labels):
    """
    Plot an iterable of slices to <path>0.png, <path>1.png etc.

    """
    pass
