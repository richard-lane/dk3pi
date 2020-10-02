"""
Script for getting real LHCb data

Creates a list of DST paths to get, and then should probably do something to get them

Run from within Ganga on lxplus

"""
import os
import re


def bk_paths(stripping_versions, mag_types, file_type) -> list:
    """
    Find the CERN grid bookkeeping paths for D->K3Pi

    Contains a lot of internal state about the allowed stripping versions, etc.

        stripping_versions: iterable of stripping versions to use as strs
        mag_types         : iterable of magnetisation directions to use as strs (i.e. "MagDown" and/or "MagUp")
        file_type         : DST file name as str, e.g. "CHARMCOMPLETEEVENT.DST"

    """
    # All possible magnet polarities
    knownMagTypes = ["MagDown", "MagUp"]

    # All possible stripping versions
    knownStrippingVersions = {
        "2011": ["Stripping21r1", "Stripping20r1"],
        "2012": ["Stripping21", "Stripping20"],
        "2015": ["Stripping24"],
        "2016": ["Stripping28"],
        "2017": ["Stripping29"],
        "2018": ["Stripping34"],
    }

    # Beam energies for each year
    beamEnergies = {
        "2011": "3500",
        "2012": "4000",
        "2015": "6500",
        "2016": "6500",
        "2017": "6500",
        "2018": "6500",
    }

    # Reco versions for each stripping version
    recoVersion = {
        "Stripping20": "14",
        "Stripping20r1": "14",
        "Stripping21": "14",
        "Stripping21r1": "14",
        "Stripping24": "15a",
        "Stripping28": "16",
        "Stripping29": "17",
        "Stripping34": "18",
    }

    # Init list of LFNs
    paths = []

    for stripping in stripping_versions:
        print("Stripping version: " + stripping)

        # Find the year corresponding to this stripping version
        # Doesn't do anything in the case where a stripping version occurs over multiple years, but hopefully this hasn't happened
        year = None
        for key in knownStrippingVersions:
            if stripping in knownStrippingVersions[key]:
                year = key
        assert year

        for magtype in mag_types:
            assert magtype in knownMagTypes
            print("\tMagnet setting: " + magtype)

            paths.append(
                "/LHCb/Collision%s/Beam%sGeV-VeloClosed-%s/Real Data/Reco%s/%s/90000000/%s"
                % (
                    year[2:],
                    beamEnergies[year],
                    magtype,
                    recoVersion[stripping],
                    stripping,
                    file_type,
                )
            )

    return paths


def check_bkfiles_exist(bookkeeping_paths: list) -> None:
    """
    Check whether an iterable of file paths all exist in the DIRAC bookkeeping

    Raises something (slowly...) if a file doesn't exist

    """
    print(
        "Checking if\n\t"
        + "\n\t".join(bookkeeping_paths)
        + "\nexist; might throw a weird error if one doesn't"
    )
    for bookkeeping_path in bookkeeping_paths:
        # Make a request that will fail if the path doesn't exist, and will succeed otherwise
        BKQuery(type="Path", dqflag="OK", path=bookkeeping_path).getDataset()


def main():
    # Choose magnet polarities + stripping versions
    # We're interested in both magnet polarities
    # We want all the stripping versions that contain charm data
    magtypes = ["MagDown", "MagUp"]
    strippings = [
        "Stripping21r1",
        "Stripping20r1",
        "Stripping21",
        "Stripping20",
        "Stripping24",
        "Stripping28",
        "Stripping29",
        "Stripping34",
    ]

    bk_files = bk_paths(strippings, magtypes, "CHARMCOMPLETEEVENT.DST")

    # Check the files exist
    check_bkfiles_exist(bk_files)

    for path in bk_files:
        # Get the data
        print("Get " + path)

# Unfortunately we can't wrap this in if name==main since we need to run it via ganga
main()


"""
A remnant from the Old Years
        bkpath = (
            "/LHCb/Collision%s/Beam%sGeV-VeloClosed-%s/Real Data/Reco%s/%s/90000000/%s"
            % (
                datatype[2:],
                beamEnergy[datatype],
                magtype,
                recoVersion[stripping],
                stripping,
                fileType,
            )
        )
        print("Using BK path: " + bkpath)
        bkq = BKQuery(type="Path", dqflag="OK", path=bkpath)
        ds = bkq.getDataset()

        # Find running period
        whichData = None
        for k, v in knownStrippingVersions.iteritems():
            if stripping in v:
                whichData = k
                break
        if not whichData:
            raise Exception("Unknown running period")

        params = {}
        params["stripping"] = stripping
        params["magtype"] = magtype
        params["whichData"] = whichData

        b = BenderModule(
            version=benderVersion, directory=userReleaseArea, module=moduleFile
        )
        b.params = params

        j = Job()
        j.name = "Collision" + datatype[2:] + "_" + magtype + "_" + stripping
        j.application = b
        j.backend = Dirac()
        j.inputdata = ds
        # j.backend.settings['BannedSites'] = ['LCG.IN2P3.fr']

        # NB remember to change the tuple name in the Bender script to match this!
        tupleFile = (
            "Bd2hhpi0-Collision"
            + datatype[2:]
            + "-"
            + magtype
            + "-"
            + stripping
            + ".root"
        )
        j.outputfiles = [DiracFile(tupleFile), LocalFile("summary.xml")]
        j.splitter = SplitByFiles(filesPerJob=50)
"""
