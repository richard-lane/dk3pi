"""
Definitions of things used to split up data; phsp bin, time bins, data taking year, magnet polarity, etc.

NB: importing this file may run attempt to build the CF and DCS wrapper libraries, if they are not already built

"""
import os
from subprocess import run, PIPE
from numpy import pi, exp

ALLOWED_MAGNET = {"MagUp", "MagDown"}
ALLOWED_YEAR = {"2016", "2017", "2018"}
DECAY_CODES = {"Phsp": 27165070, "RS": 27165071, "WS": 27165072}

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
DCS_OFFSET = OFFSET_MAG * exp((0 + 1j) * OFFSET_PHASE * pi / 180.0)

MODEL_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "amplitude_models")
CF_LIB = os.path.abspath(os.path.join(MODEL_DIR, "cf_wrapper.so"))
DCS_LIB = os.path.abspath(os.path.join(MODEL_DIR, "dcs_wrapper.so"))

BDT_PATHS = os.path.abspath(os.path.join(os.path.dirname(__file__), "bdt_paths.pickle"))

# Tree names for WG productions
RS_TREE = "Hlt2Dstp2D0Pip_D02KmPimPipPip_Tuple/DecayTree"
WS_TREE = "Hlt2Dstp2D0Pip_D02KpPimPimPip_Tuple/DecayTree"

# Paths to data generated with models, stored on lxplus.cern.ch
RS_AMPGEN_PATH = "/eos/home-r/rilane/Documents/data/ampgen_Dbar_RS.root"
WS_AMPGEN_PATH = "/eos/home-r/rilane/Documents/data/ampgen_WS.root"

# The function handles used to evaluate amplitudes
# Will be set at run-time and cached here
DCS_FCN = None
CF_FCN = None


def build_libs():
    """
    Build the AmpGen wrapper libraries

    """
    build_script = os.path.join(MODEL_DIR, "build.sh")
    print(
        f"Building AmpGen wrapper libs, required from \n\t{__file__} ...",
        end="",
        flush=True,
    )
    run([build_script], stdout=PIPE, stderr=PIPE, check=True)
    print("done")


# If our shared libraries haven't been built, build them
if not os.path.exists(CF_LIB) or not os.path.exists(DCS_LIB):
    build_libs()


def assign_dcs(handle):
    global DCS_FCN
    DCS_FCN = handle


def assign_cf(handle):
    global CF_FCN
    CF_FCN = handle
