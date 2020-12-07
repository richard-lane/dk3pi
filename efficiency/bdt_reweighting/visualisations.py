"""
Helpers for visualising our multidimensional phase space distributions

"""
import matplotlib.pyplot as plt
import numpy as np


class Slices:
    def __init__(self, num_slices, num_bins, bin_limits, plot_index, slice_index):
        """
        Init an instance for plotting histogram projections of the plot_index'th parameter,
        taking slices along the slice_index'th parameter

        bin_limits should be an iterable (lower_limit, upper_limit)

        """
        self._plot_index = plot_index
        self._slice_index = slice_index

        # Create an array for the limits of our histogram slices
        self._bin_limits = np.linspace(bin_limits[0], bin_limits[1], num_bins + 1)

        # Create an array for the limits of our slices
        self._slice_limits = np.linspace(bin_limits[0], bin_limits[1], num_slices + 1)

        # Create an array of arrays for the bin contents of each slice
        self._slice_array = np.zeros((num_slices, num_bins), dtype=np.float64)

    def add_point(self, point, weight=1):
        """
        Add a point to the correct bin in the correct histogram

        """
        # Find the right histogram
        # If the point's _slice_index'th co-ordinate is in the N'th bin of _slice_limits,
        # it belongs in the N'th histogram slice
        this_slice = np.digitize(point[self._slice_index], self._slice_limits) - 1
        assert 0 <= this_slice < len(self._slice_limits)

        # Find the right bin in this histogram
        this_bin = np.digitize(point[self._plot_index], self._bin_limits) - 1
        assert 0 <= this_bin < len(self._bin_limits)

        # Increment the right bin
        self._slice_array[this_slice][this_bin] += weight
        self._num_points += weight

    def add_points(self, points, weights=None):
        """
        Add a series of points to the correct bins in the correct histograms

        """
        w = weights if weights is not None else np.ones_like(points[:, 0])
        assert len(points) == len(w)

        for point, weight in zip(points, weights):
            self.add_point(point, weight)

    # Array of histograms representing each slice through phase space
    _slice_array = []

    # Bin limits for finding which bin a point belongs in
    _bin_limits = []

    # Bin limits for finding which histogram a point belongs in
    _slice_limits = []

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
    assert len(slices) == len(labels)

    num_slices = len(slices[0]._slice_array)
    slice_limits = slices[0]._slice_limits
    bin_limits = slices[0]._bin_limits
    for s in slices:
        assert len(s._slice_array) == num_slices
        assert np.allclose(s._slice_limits, slice_limits)
        assert np.allclose(s._bin_limits, bin_limits)

    # Find the bin centres for plotting
    centres = np.mean(np.vstack([bin_limits[0:-1], bin_limits[1:]]), axis=0)

    for i in range(num_slices):
        slice_path = path + str(i) + ".png"

        for s, label in zip(slices, labels):
            # Scale such that the sum of all slices have an area of 1
            counts = np.copy(s._slice_array[i]) / s._num_points

            # Plot
            plt.plot(centres, counts, label=label)
        plt.legend()
        plt.savefig(slice_path)
        plt.clf()
