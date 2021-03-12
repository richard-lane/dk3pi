"""
General utilities

"""


def _unsafe_bin(data, bins):
    """
    Bin a data point into bins without bounds checking.
    Will return 0 if the point underflows; None if overflow

    :param data: data point
    :param bins: left edge of each bin, and rightmost edge of highest bin
    :returns: bin number, from 0

    """
    for i, edge in enumerate(bins[1:]):
        if data < edge:
            return i


def _bin(data, bins):
    """
    Bin a data point into bins

    :param data: data point
    :param bins: left edge of each bin, and rightmost edge of highest bin
    :returns: bin number, from 0
    :raises ValueError: if data out of range

    """
    if data < bins[0] or data > bins[-1]:
        raise ValueError(
            f"{data} out of range; bins cover range [{bins[0]}, {bins[-1]}]"
        )

    return _unsafe_bin(data, bins)
