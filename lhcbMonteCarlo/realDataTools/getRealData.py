"""
Get the real data

"""
import os
import re


def bk_paths(stripping_versions, mag_types) -> list:
    """
    Find the CERN grid bookkeeping paths for D->K3Pi

    Takes in an iterable of stripping versions and magnet polarities to use; returns a list of LFNs

    Contains a lot of internal state about the allowed stripping versions, etc.

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

        # Find the year corresponding to this stripping version
        # Doesn't do anything in the case where a stripping version occurs over multiple years, but hopefully this hasn't happened
        year = None
        for key in knownStrippingVersions:
            if stripping in knownStrippingVersions[key]:
                year = key
        assert year

        print("Creating job(s) for stripping: " + stripping)

        fileType = "CHARMCOMPLETEEVENT.DST"
        for magtype in magtypes:
            assert magtype in knownMagTypes
            print("With magnet setting: " + magtype)

            paths.append(
                "/LHCb/Collision%s/Beam%sGeV-VeloClosed-%s/Real Data/Reco%s/%s/90000000/%s"
                % (
                    year[2:],
                    beamEnergies[year],
                    magtype,
                    recoVersion[stripping],
                    stripping,
                    fileType,
                )
            )

    return paths


def check_bkfile_exists(bookkeeping_path: str) -> None:
    """
    Check whether a file exists in the DIRAC bookkeeping

    Raises something (slowly...) if the file doesn't exist

    """
    print(
        "Checking if "
        + bookkeeping_path
        + " exists; might throw a weird error if it doesn't"
    )
    BKQuery(type="Path", dqflag="OK", path=bookkeeping_path).getDataset()


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

bk_files = bk_paths(strippings, magtypes)

for path in bk_files:
    # Check the file exists
    check_bkfile_exists(path)

    # Get the data
    print("Get " + path)

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
