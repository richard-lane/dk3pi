"""
Python binding for Z scan with CLEO fitter

"""
import numpy as np
import sys, os

def _mungepath():
    """
    Add the right location of our CLEO C++ lib to python path
    """
    from pathlib import PurePath

    p = PurePath(os.path.dirname(os.path.abspath(__file__)))

    new_path = list(p.parts)
    new_path.insert(-2, "build")

    sys.path.append(str(PurePath("").joinpath(*new_path)))
_mungepath()

import libcleoScan

def cleoLikelihoods(reZVals: np.ndarray,
             imZVals: np.ndarray,
             decayParams: np.ndarray,
             binNumber: int) -> np.ndarray:
    """
    Evaluate CLEO likelihoods, fixing Re(Z) and Im(Z) to the provided values; return an array of likelihoods

    Return value indexed as (i + r * nImZVals) for imaginary+real (i, r)

    """
    return np.array(libcleoScan.cleoLikelihoods(reZVals,
                                      imZVals,
                                      decayParams,
                                      binNumber))

def charmScan(reZVals: np.ndarray,
              imZVals: np.ndarray,
              rsDecayTimes: np.ndarray,
              rsWeights: np.ndarray,
              wsDecayTimes: np.ndarray,
              wsWeights: np.ndarray,
              binLimits: np.ndarray,
              initialVals: np.ndarray,
              initialErrs: np.ndarray) -> np.ndarray:
    """
    Python wrapper around my C++ lib so I can tell what's going on

    Perform fits fixing Re(Z) and Im(Z) to the provided values; return an array of likelihoods

    Return value indexed as (i + r * nImZVals) for imaginary+real (i, r)

    """
    return np.array(libcleoScan.charmLikelihoods(reZVals,
                                      imZVals,
                                      rsDecayTimes,
                                      rsWeights,
                                      wsDecayTimes,
                                      wsWeights,
                                      binLimits,
                                      initialVals,
                                      initialErrs))

def cleoScan(reZVals: np.ndarray,
             imZVals: np.ndarray,
             rsDecayTimes: np.ndarray,
             rsWeights: np.ndarray,
             wsDecayTimes: np.ndarray,
             wsWeights: np.ndarray,
             binLimits: np.ndarray,
             initialVals: np.ndarray,
             initialErrs: np.ndarray,
             binNumber: int) -> np.ndarray:
    """
    Python wrapper around my C++ lib so I can tell what's going on

    Perform fits fixing Re(Z) and Im(Z) to the provided values; return an array of likelihoods
    
    Return value indexed as (i + r * nImZVals) for imaginary+real (i, r)

    """
    return np.array(libcleoScan.cleoZScan(reZVals,
                                      imZVals,
                                      rsDecayTimes,
                                      rsWeights,
                                      wsDecayTimes,
                                      wsWeights,
                                      binLimits,
                                      initialVals,
                                      initialErrs,
                                      binNumber))

