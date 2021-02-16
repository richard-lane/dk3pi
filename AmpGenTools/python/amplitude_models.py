"""
TODO make these better, we shouldn't have to read from the shared lib every time we call the function

"""
import os
import subprocess
import ctypes
import numpy as np


# If our shared libraries haven't been built, build them
if not os.path.exists(
    os.path.abspath(os.path.join(os.path.dirname(__file__), "cf_wrapper.so"))
) or not os.path.exists(
    os.path.abspath(os.path.join(os.path.dirname(__file__), "dcs_wrapper.so"))
):
    build_script = os.path.join(os.path.dirname(__file__), "build.sh")
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


def cf_amplitude(event: np.ndarray, k_charge: int):
    """
    Find the CF amplitude of an event

    :param event   : a 16-element array of kinematic parameters (kpx, kpy, kpz, kE, ...) for k, pi1, pi2, pi3
    :param k_charge: the charge of the kaon, i.e. +1 or -1
    :raises        : many things
    :returns       : the complex amplitude as an instance of the builtin complex

    """
    # Find the function that we need in our shared library
    lib = os.path.abspath(os.path.join(os.path.dirname(__file__), "cf_wrapper.so"))
    cf_fcn = _find_fcn(lib, "cf_wrapper")
    # Need to copy the event otherwise our array is not contiguous, and so can't be read by ctypes
    event = np.copy(event)

    # Convert our arguments to compatible types
    event_p = event.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    k_charge_ctype = ctypes.c_int(k_charge)

    # Call the function
    retval = cf_fcn(event_p, k_charge_ctype)

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
    lib = os.path.abspath(os.path.join(os.path.dirname(__file__), "dcs_wrapper.so"))
    dcs_fcn = _find_fcn(lib, "dcs_wrapper")
    # Need to copy the event otherwise our array is not contiguous, and so can't be read by ctypes
    event = np.copy(event)

    # Convert our arguments to compatible types
    event_p = event.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    k_charge_ctype = ctypes.c_int(k_charge)

    # Call the function
    retval = dcs_fcn(event_p, k_charge_ctype)

    return complex(retval.real, retval.imag)
