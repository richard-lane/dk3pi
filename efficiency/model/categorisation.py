"""
Definitions of things used to split up data; phsp bin, time bins, data taking year, magnet polarity, etc.

NB: importing this file may run attempt to build the CF and DCS wrapper libraries, if they are not already built

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
D_LIFETIME = 0.00041  # nanoseconds
TIME_BINS = (-1.0, 0.0, 0.94, 1.185, 1.40, 1.62, 1.85, 2.13, 2.45, 2.87, 3.5, 8.0, 19.0)

# We will veto any events with M(pipi) in this range, for any pair of pions
KS_MASS = 0.497614  # MeV
VETO_WIDTH = 0.010

# If our shared libraries haven't been built, build them
model_dir = os.path.join(os.path.dirname(__file__), "amplitude_models")
if not os.path.exists(
    os.path.abspath(os.path.join(model_dir, "cf_wrapper.so"))
) or not os.path.exists(os.path.abspath(os.path.join(model_dir, "dcs_wrapper.so"))):
    build_script = os.path.join(model_dir, "build.sh")
    print(
        f"Building AmpGen wrapper libs, required from \n\t{__file__} ...",
        end="",
        flush=True,
    )
    subprocess.run(
        [build_script], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
    )
    print("done")


def time_bin(time):
    """
    Find which time bin an event belongs in

    :param time: decay time in nanoseconds
    :returns: Index of the required time bin

    """
    lifetimes = time / D_LIFETIME

    # Over/underflow
    if lifetimes < TIME_BINS[0] or lifetimes > TIME_BINS[-1]:
        raise ValueError(
            f"Time {time}ns out of range; time bins cover range [{TIME_BINS[0]}, {TIME_BINS[-1]}] lifetimes"
        )

    for i, edge in enumerate(TIME_BINS[1:]):
        if lifetimes < edge:
            return i


def phsp_bin(event):
    """
    Find which phase space bin an event belongs in

    """
    ...


def vetoed(mass):
    """
    Check whether a mass is within the veto range of the Ks

    Intended to check the mass of a (pipi) system

    :param: invariant mass in MeV
    :returns: bool; whether the event is veto'd

    """
    return abs(mass - KS_MASS) < VETO_WIDTH
