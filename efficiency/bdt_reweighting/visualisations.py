"""
Helpers for visualising our multidimensional phase space distributions

"""


class Slices:
    def __init__(self, num_slices, num_bins, bin_limits, plot_index, slice_index):
        """
        Init an instance for plotting histogram projections of the plot_index'th parameter,
        taking slices along the slice_index'th parameter

        bin_limits should be an iterable (lower_limit, upper_limit)

        """
        _plot_index = plot_index
        _slice_index = slice_index

        # Create an array for the limits of our histogram slices
        _bin_limits = np.linspace(bin_limits[0], bin_limits[1], num_bins + 1)

        # Create an array for the limits of our slices
        _slice_limits = np.linspace(bin_limits[0], bin_limits[1], num_slices + 1)

        # Create an array of arrays for the bin contents of each slice
        _slice_array = np.zeros((num_slices, num_bins))

    def add_point(self, point, weight=1):
        """
        Add a point to the correct bin in the correct histogram

        """
        pass

    def add_points(self, points, weights=None):
        """
        Add a series of points to the correct bins in the correct histograms

        """
        pass

    # Array of histograms representing each slice through phase space
    _slice_array = None

    # Bin limits for finding which bin a point belongs in
    _bin_limits = None

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
