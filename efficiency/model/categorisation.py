from . import definitions

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
    return _bin(phase, definitions.PHSP_BINS)


def cf_amplitude(event: np.ndarray, k_charge: int):
    """
    Find the CF amplitude of an event

    :param event   : a 16-element array of kinematic parameters (kpx, kpy, kpz, kE, ...) for k, pi1, pi2, pi3
    :param k_charge: the charge of the kaon, i.e. +1 or -1
    :raises        : many things
    :returns       : the complex amplitude as an instance of the builtin complex

    """
    # Find the function that we need in our shared library
    cf_fcn = _find_fcn(defnitionsCF_LIB, "cf_wrapper")

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
    dcs_fcn = _find_fcn(definitions.DCS_LIB)

    retval = _amplitude(dcs_fcn, event, k_charge)
    return complex(retval.real, retval.imag)


def relative_phase(event: np.ndarray, k_charge: int):
    """
    Find the relative phase between the DCS and CF models for a given event, subject to the dcs model offset
    :param event   : a 16-element array of kinematic parameters (kpx, kpy, kpz, kE, ...) for k, pi1, pi2, pi3
    :param k_charge: the charge of the kaon, i.e. +1 or -1
    :returns       : relative phase in degrees

    """
    cf = cf_amplitude(event, k_charge)
    dcs = dcs_amplitude(event, k_charge) * definitions.DCS_OFFSET

    return np.angle(cf * dcs.conjugate(), deg=True)


def time_bin(time):
    """
    Find which time bin an event belongs in

    :param time: decay time in nanoseconds
    :returns: Index of the required time bin

    """
    lifetimes = time / definitions.D_LIFETIME

    return _bin(lifetimes, definitions.TIME_BINS)


def phsp_bin(event, k_charge):
    """
    Find which phase space bin an event belongs in

    """
    return _bin(relative_phase(event, k_charge), definitions.PHSP_BINS)


def vetoed(mass):
    """
    Check whether a mass is within the veto range of the Ks

    Intended to check the mass of a (pipi) system

    :param: invariant mass in MeV
    :returns: bool; whether the event is veto'd

    """
    return abs(mass - definitions.KS_MASS) < definitions.VETO_WIDTH
