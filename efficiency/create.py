from model.creation import create_bdts
from model.creation import productions

from model.util import definitions

import glob
import uproot
import numpy as np


def main():
    # Read RS 2018 MagUp MC
    (
        kpx,
        kpy,
        kpz,
        ke,
        pi1px,
        pi1py,
        pi1pz,
        pi1e,
        pi2px,
        pi2py,
        pi2pz,
        pi2e,
        pi3px,
        pi3py,
        pi3pz,
        pi3e,
    ) = ([] for _ in range(16))

    # TODO SPLIT BY K CHARGES

    for filename in glob.glob(productions.get("MagDown", 2018, "RS")):
        tree = uproot.open(filename)[definitions.RS_TREE]
        print(filename)

        "D0_TAU"

        "D0_P0_PX"
        "D0_P0_PY"
        "D0_P0_PZ"
        "D0_P0_E"

        "D0_P1_PX"
        "D0_P1_PY"
        "D0_P1_PZ"
        "D0_P1_E"
        "D0_P2_PX"
        "D0_P2_PY"
        "D0_P2_PZ"
        "D0_P2_E"
        "D0_P3_PX"
        "D0_P3_PY"
        "D0_P3_PZ"
        "D0_P3_E"

    # Categorise events into phsp bins

    # Read AmpGen
    "_1_K~_Px"
    "_1_K~_Py"
    "_1_K~_Pz"
    "_1_K~_E"
    "_2_pi#_Px"
    "_2_pi#_Py"
    "_2_pi#_Pz"
    "_2_pi#_E"
    "_3_pi#_Px"
    "_3_pi#_Py"
    "_3_pi#_Pz"
    "_3_pi#_E"
    "_4_pi~_Px"
    "_4_pi~_Py"
    "_4_pi~_Pz"
    "_4_pi~_E"


if __name__ == "__main__":
    main()
