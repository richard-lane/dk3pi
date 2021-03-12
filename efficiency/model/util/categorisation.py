from . import definitions
import ctypes
import numpy as np


def time_bin(time):
    """
    Find which time bin an event belongs in

    :param time: decay time in nanoseconds
    :returns: Index of the required time bin

    """
    lifetimes = time / definitions.D_LIFETIME

    return _bin(lifetimes, definitions.TIME_BINS)


def vetoed(mass):
    """
    Check whether a mass is within the veto range of the Ks

    Intended to check the mass of a (pipi) system

    :param: invariant mass in MeV
    :returns: bool; whether the event is veto'd

    """
    return abs(mass - definitions.KS_MASS) < definitions.VETO_WIDTH
