"""
Script for getting real LHCb data

Creates a list of DST paths to get, and then should probably do something to get them

Run from within Ganga on lxplus

"""
import os
import re


def bk_paths(years, mag_types, dst_name) -> list:
    """
    Find the CERN grid bookkeeping paths for D->K3Pi

    Contains a lot of hard-coded internal state about the allowed stripping versions, etc. that I just got by looking on the LHCb stripping project website

        years    : iterable of data-taking years to consider, as strs
        mag_types: iterable of magnetisation directions to use as strs (i.e. "MagDown" and/or "MagUp")
        dst_name : DST file name as str, e.g. "CHARMCOMPLETEEVENT.DST"

    """
    # One stripping line per year that contains the data I want
    stripping_lines = {
        "2011": "StrippingDstarPromptWithD02HHHHLine",
        "2012": "StrippingDstarPromptWithD02HHHHLine",
        "2015": "StrippingDstarD2HHHHDstarD2KPiPiPiLine ",
        "2016": "StrippingDstarD2HHHHDstarD2KPiPiPiLine ",
        "2017": "StrippingDstarD2HHHHDstarD2KPiPiPiLine",
        "2018": "StrippingDstarD2HHHHDstarD2KPiPiPiLine",
    }

    # Stripping versions that contain the stripping lines I want
    # I found these by going on the LHCb stripping project site
    stripping_versions = {
        "2011": "Stripping21r1",
        "2012": "Stripping21",
        "2015": "Stripping24r2",
        "2016": "Stripping28r2",
        "2017": "Stripping29r2",
        "2018": "Stripping34",
    }

    # Beam energies
    beam_energies = {
        "2011": "3500",
        "2012": "4000",
        "2015": "6500",
        "2016": "6500",
        "2017": "6500",
        "2018": "6500",
    }

    # Reco versions for each stripping version
    # Could just as well define these per year, but this is perhaps more intuitive as a reco version is tied to a stripping version
    # I got these from https://twiki.cern.ch/twiki/bin/view/Main/ProcessingPasses, which may be obsolete if you are reading this in the future
    # If it is you can find the reco versions from the DIRAC bookkeeping browser
    reco_version = {
        "Stripping21r1": "14",
        "Stripping21": "14",
        "Stripping24r2": "15a",
        "Stripping28r2": "16",
        "Stripping29r2": "17",
        "Stripping34": "18",
    }

    paths = []
    for year in years:
        print("Year: " + year)

        # Find the year corresponding to this stripping version
        # Doesn't do anything in the case where a stripping version occurs over multiple years, but hopefully this hasn't happened
        # year = None
        # for key in knownStrippingVersions:
        #    if stripping in knownStrippingVersions[key]:
        #        year = key
        # assert year

        stripping_version = stripping_versions[year]
        for magtype in mag_types:
            print("\tMagnet setting: " + magtype)

            paths.append(
                "/LHCb/Collision%s/Beam%sGeV-VeloClosed-%s/Real Data/Reco%s/%s/90000000/%s"
                % (
                    year[2:],
                    beam_energies[year],
                    magtype,
                    reco_version[stripping_version],
                    stripping_version,
                    dst_name,
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
    # We're interested in both magnet polarities
    magtypes = ("MagDown", "MagUp")

    # We want data from all years
    years = ("2011", "2012", "2015", "2016", "2017", "2018")

    bk_files = bk_paths(years, magtypes, "CHARMCOMPLETEEVENT.DST")

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
