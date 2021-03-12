from . import create_bdts
from . import definitions
from . import productions
import glob
import uproot


def main():
    # Read RS 2018 MagUp MC
    kpx, kpy, kpz, ke, pi1px, pi1py, pi1pz, pi1e, pi2px, pi2py, pi2pz, pi2e, pi3px, pi3py, pi3pz, pi3e = (
        [] for _ in range(16)
    )

    # TODO SPLIT BY K CHARGES

    for filename in glob.glob(productions.production[2018]["MagUp"][27165071]):
        tree = uproot.open(filename)[definitions.RS_TREE]

        t.extend(tree.array("D0_TAU"))

        kpx.extend(tree.array("D0_P0_PX"))
        kpy.extend(tree.array("D0_P0_PY"))
        kpz.extend(tree.array("D0_P0_PZ"))
        ke.extend(tree.array("D0_P0_E"))

        pi1px.extend(tree.array("D0_P1_PX"))
        pi1px.extend(tree.array("D0_P1_PY"))
        pi1px.extend(tree.array("D0_P1_PZ"))
        pi1e.extend(tree.array("D0_P1_E"))
        pi2px.extend(tree.array("D0_P2_PX"))
        pi2px.extend(tree.array("D0_P2_PY"))
        pi2px.extend(tree.array("D0_P2_PZ"))
        pi2e.extend(tree.array("D0_P2_E"))
        pi3px.extend(tree.array("D0_P3_PX"))
        pi3px.extend(tree.array("D0_P3_PY"))
        pi3px.extend(tree.array("D0_P3_PZ"))
        pi3e.extend(tree.array("D0_P3_E"))

    # Categorise events into phsp bins
    bin_numbers = -1 * np.ones(len(t))
    for i in range(len(t)):
        evt = np.array(
            (
                kpx[i],
                kpy[i],
                kpz[i],
                ke[i],
                pi1px[i],
                pi1py[i],
                pi1pz[i],
                pi1e[i],
                pi2px[i],
                pi2py[i],
                pi2pz[i],
                pi2e[i],
                pi3px[i],
                pi3py[i],
                pi3pz[i],
                pi3e[i],
            )
        )
        bin_numbers[i] = categorisation.phsp_bin(evt, +1)

    # Read AmpGen
    tree = uproot.open(filename)[definitions.RS_TREE]
    ag_kpx, ag_kpy, ag_kpz, ag_ke, ag_pi1px, ag_pi1py, ag_pi1pz, ag_pi1e, ag_pi2px, ag_pi2py, ag_pi2pz, ag_pi2e, ag_pi3px, ag_pi3py, ag_pi3pz, ag_pi3e = (
        tree.array(branch)
        for branch in (
            "_1_K~_Px",
            "_1_K~_Py",
            "_1_K~_Pz",
            "_1_K~_E",
            "_2_pi#_Px",
            "_2_pi#_Py",
            "_2_pi#_Pz",
            "_2_pi#_E",
            "_3_pi#_Px",
            "_3_pi#_Py",
            "_3_pi#_Pz",
            "_3_pi#_E",
            "_4_pi~_Px",
            "_4_pi~_Py",
            "_4_pi~_Pz",
            "_4_pi~_E",
        )
    )

    ag_t.extend(tree.array("D0_TAU"))

    ag_kpx.extend(tree.array("D0_P0_PX"))
    ag_kpy.extend(tree.array("D0_P0_PY"))
    ag_kpz.extend(tree.array("D0_P0_PZ"))
    ag_ke.extend(tree.array("D0_P0_E"))

    ag_pi1px.extend(tree.array("D0_P1_PX"))
    ag_pi1px.extend(tree.array("D0_P1_PY"))
    ag_pi1px.extend(tree.array("D0_P1_PZ"))
    ag_pi1e.extend(tree.array("D0_P1_E"))
    ag_pi2px.extend(tree.array("D0_P2_PX"))
    ag_pi2px.extend(tree.array("D0_P2_PY"))
    ag_pi2px.extend(tree.array("D0_P2_PZ"))
    ag_pi2e.extend(tree.array("D0_P2_E"))
    ag_pi3px.extend(tree.array("D0_P3_PX"))
    ag_pi3px.extend(tree.array("D0_P3_PY"))
    ag_pi3px.extend(tree.array("D0_P3_PZ"))
    ag_pi3e.extend(tree.array("D0_P3_E"))

    # Categorise events into phsp bins
    ag_bin_numbers = -1 * np.ones(len(t))
    for i in range(len(t)):
        evt = np.array(
            (
                ag_kpx[i],
                ag_kpy[i],
                ag_kpz[i],
                ag_ke[i],
                ag_pi1px[i],
                ag_pi1py[i],
                ag_pi1pz[i],
                ag_pi1e[i],
                ag_pi2px[i],
                ag_pi2py[i],
                ag_pi2pz[i],
                ag_pi2e[i],
                ag_pi3px[i],
                ag_pi3py[i],
                ag_pi3pz[i],
                ag_pi3e[i],
            )
        )
        ag_bin_numbers[i] = categorisation.phsp_bin(evt, +1)

    # For each phsp bin, train a BDT and pickle it at a location
    for i, path in enumerate(("path1", "path2", "path3", "path4")):
        target = 0
        original = 0

        create_bdts.train(
            target, original, path, 2018, "MagUp", i, definitions.DECAY_CODES["RS"]
        )


if __name__ == "__main__":
    main()
