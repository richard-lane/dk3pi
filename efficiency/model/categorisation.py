"""
Definitions of things used to split up data; phsp bin, time bins, data taking year, magnet polarity, etc.

NB: importing this file may run attempt to build the CF and DCS wrapper libraries, if they are not already built

TODO make these better, we shouldn't have to read from the shared lib every time we call the function

"""
import os
import subprocess
import ctypes
import numpy as np

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

# Scale + rotate amplitudes so that dcs/cf amplitude ratio ~ 0.055 and relative strong phase ~ 0
OFFSET_MAG = 0.0601387
OFFSET_PHASE = 1.04827  # degrees
DCS_OFFSET = OFFSET_MAG * np.exp((0 + 1j) * OFFSET_PHASE * np.pi / 180.0)

MODEL_DIR = os.path.join(os.path.dirname(__file__), "amplitude_models")
CF_LIB = os.path.abspath(os.path.join(MODEL_DIR, "cf_wrapper.so"))
DCS_LIB = os.path.abspath(os.path.join(MODEL_DIR, "dcs_wrapper.so"))

# If our shared libraries haven't been built, build them
if not os.path.exists(CF_LIB) or not os.path.exists(DCS_LIB):
    build_script = os.path.join(MODEL_DIR, "build.sh")
    print(
        f"Building AmpGen wrapper libs, required from \n\t{__file__} ...",
        end="",
        flush=True,
    )
    subprocess.run(
        [build_script], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
    )
    print("done")


class Complex(ctypes.Structure):
    """
    Class that enables us to read a simple struct returned from C

    """

    _fields_ = [("real", ctypes.c_double), ("imag", ctypes.c_double)]


def _find_fcn(lib_name, fcn_name):
    """
    Find the named function in a shared library

    This will look in the shared library for the function of the provided name, and enforce that its
    argument types are an array and an integer, and that it returns a Complex

    :param lib_name: Location of the shared library to look in
    :param fcn_name: Name of the function to return
    :return        : A ctypes handle to the function

    """
    dll = ctypes.cdll.LoadLibrary(lib_name)
    fcn = getattr(dll, fcn_name)
    fcn.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int]
    fcn.restype = Complex

    return fcn


def _amplitude(fcn, event, k_charge):
    """
    Find an amplitude from an event and amplitude fcn

    """
    # Need to copy the event otherwise our array may not be contiguous, and so can't be read by ctypes
    event = np.copy(event)

    # Convert our arguments to compatible types
    event_p = event.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    k_charge_ctype = ctypes.c_int(k_charge)

    # Call the function
    return fcn(event_p, k_charge_ctype)


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

    for i, edge in enumerate(bins[1:]):
        if data < edge:
            return i


def _bin_from_phase(phase):
    """
    Bin a phase in degrees into PHSP_BINS

    Returns bin numbered from 0

    """
    return _bin(phase, PHSP_BINS)


def cf_amplitude(event: np.ndarray, k_charge: int):
    """
    Find the CF amplitude of an event

    :param event   : a 16-element array of kinematic parameters (kpx, kpy, kpz, kE, ...) for k, pi1, pi2, pi3
    :param k_charge: the charge of the kaon, i.e. +1 or -1
    :raises        : many things
    :returns       : the complex amplitude as an instance of the builtin complex

    """
    # Find the function that we need in our shared library
    cf_fcn = _find_fcn(CF_LIB, "cf_wrapper")

    retval = _amplitude(cf_fcn, event, k_charge)
    return complex(retval.real, retval.imag)


def dcs_amplitude(event: np.ndarray, k_charge: int):
    """
    Find the DCS amplitude of an event

    :param event   : a 16-element array of kinematic parameters (kpx, kpy, kpz, kE, ...) for k, pi1, pi2, pi3
    :param k_charge: the charge of the kaon, i.e. +1 or -1
    :raises        : many things
    :returns       : the complex amplitude as an instance of the builtin complex

    """
    # Find the function that we need in our shared library
    dcs_fcn = _find_fcn(DCS_LIB)

    retval = _amplitude(dcs_fcn, event, k_charge)
    return complex(retval.real, retval.imag)


def time_bin(time):
    """
    Find which time bin an event belongs in

    :param time: decay time in nanoseconds
    :returns: Index of the required time bin

    """
    lifetimes = time / D_LIFETIME

    return _bin(lifetimes, TIME_BINS)


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
