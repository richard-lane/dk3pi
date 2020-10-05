"""
Script for getting real LHCb data

Creates a list of DST paths to get, and then should probably do something to get them

Run from within Ganga on lxplus

"""
import os
import re


def bk_paths(years, mag_types, dst_name) -> dict:
    """
    Find the CERN grid bookkeeping paths for D->K3Pi

    Contains a lot of hard-coded internal state about the allowed stripping versions, etc. that I just got by looking on the LHCb stripping project website

        years    : iterable of data-taking years to consider, as strs
        mag_types: iterable of magnetisation directions to use as strs (i.e. "MagDown" and/or "MagUp")
        dst_name : DST file name as str, e.g. "CHARMCOMPLETEEVENT.DST"

    Returns a dict of {year: [DST paths]}
    Years are strs

    """
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

    paths = dict()
    for year in years:
        paths[year] = []
        print("Year: " + year)

        stripping_version = stripping_versions[year]
        for magtype in mag_types:
            print("\tMagnet setting: " + magtype)
            paths[year].append(
                f"/LHCb/Collision{year[2:]}/Beam{beam_energies[year]}GeV-VeloClosed-{magtype}/"
                f"Real Data/Reco{reco_version[stripping_version]}/{stripping_version}/90000000/{dst_name}"
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
        print(bookkeeping_path + " exists")


def submit_job(dst_paths: list, stripping_lines) -> None:
    """
    Submit a job to the grid; config defined in ./nTupleOptions.py

        dst_path      : list of DST locations in the LHCb bookkeeping
        stripping_line: iterable of str reprs of the LHCb stripping line to use for this dataset

    """
    # Check that the DSTs we're looking for all exist
    check_bkfiles_exist(dst_paths)

    # If the directory where we will be storing my local copy of DaVinci already exists, raise
    if os.path.exists("./DaVinciDev_v45r1/"):
        raise Exception("rm the davnci dir before submitting job")

    # Init DaVinci
    myApp = prepareGaudiExec("DaVinci", "v45r1", myPath=".")

    # Init job + tell Ganga where my config file is
    j = Job(name="Prompt DK3Pi Data")
    j.application = myApp
    j.application.options = ["nTupleOptions.py"]
    j.application.platform = "x86_64-centos7-gcc8-opt"

    # Tell Ganga what data to use
    assert (
        False
    ), "idk how to tell DaVinci about all our files/stripping lines in the right way"
    j.inputdata = dst_paths

    # Job uses Dirac backend i guess (?)
    j.backend = Dirac()

    # Tell the job to split itself up a bit
    j.splitter = SplitByFiles(filesPerJob=50)

    # Configure where to send the output
    j.outputfiles = [DiracFile("*.root")]

    # Submit job
    j.submit()


def main():
    # We want data from all years
    years = ("2011", "2012", "2015", "2016", "2017", "2018")
    bookkeeping_paths = bk_paths(years, ("MagDown", "MagUp"), "CHARMCOMPLETEEVENT.DST")

    # One stripping line per year that contains the data I want
    # From the LHCb stripping project site
    stripping_lines = {
        "2011": "StrippingDstarPromptWithD02HHHHLine",
        "2012": "StrippingDstarPromptWithD02HHHHLine",
        "2015": "StrippingDstarD2HHHHDstarD2KPiPiPiLine ",
        "2016": "StrippingDstarD2HHHHDstarD2KPiPiPiLine ",
        "2017": "StrippingDstarD2HHHHDstarD2KPiPiPiLine",
        "2018": "StrippingDstarD2HHHHDstarD2KPiPiPiLine",
    }

    for key in bookkeeping_paths:
        # Create a job for one bk path to test
        year = "2011"
        print(bookkeeping_paths[year])
        submit_job(bookkeeping_paths[year][0:1], stripping_lines[year])


# Unfortunately we can't wrap this in if name==main since we need to run it via ganga
main()

