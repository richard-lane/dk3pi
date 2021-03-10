"""
Definitions of things used to split up data; phsp bin, time bins, data taking year, magnet polarity, etc.

"""
ALLOWED_MAGNET = {"MagUp", "MagDown"}
ALLOWED_YEAR = {
    "2011",
    "2012",
    "2013",
    "2014",
    "2015",
    "2016",
    "2017",
    "2018",
}  # Check these

# Phsp bin boundaries in degrees; defined using Tim's amplitude models
PHSP_BINS = (-180.0, -39.0, 0.0, 43.0, 180.0)

# Time bin boundaries in D-lifetimes
TIME_BINS = (-1.0, 0.0, 0.94, 1.185, 1.40, 1.62, 1.85, 2.13, 2.45, 2.87, 3.5, 8.0, 19.0)

# Veto any events with M(pipi) in this range, for any pair of pions
KS_VETO = (0.487614, 0.507614)
