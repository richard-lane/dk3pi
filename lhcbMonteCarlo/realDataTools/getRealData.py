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
        stripping_version = stripping_versions[year]
        for magtype in mag_types:
            paths[year].append(
                f"/LHCb/Collision{year[2:]}/Beam{beam_energies[year]}GeV-VeloClosed-{magtype}/"
                f"Real Data/Reco{reco_version[stripping_version]}/{stripping_version}/90000000/{dst_name}"
            )

    return paths


def check_bkfile_exists(bookkeeping_path: str) -> None:
    """
    Check whether a file path exists in the DIRAC bookkeeping

    Raises something (slowly...) if a file doesn't exist

    """
    # Make a request that will fail if the path doesn't exist, and will succeed otherwise
    BKQuery(type="Path", dqflag="OK", path=bookkeeping_path).getDataset()
    print(bookkeeping_path + " exists")


def create_davinci_application(path: str, version: str):
    """
    Create a local copy of the Davinic application

    returns whatever prepareGaudiExec() returns, idk i can't find any documentation on it ???

    """
    if os.path.exists(os.path.join(path, "DaVinciDev_v45r1/")):
        raise Exception("rm the davnci dir before submitting job")

    return prepareGaudiExec("DaVinci", version, myPath=path)


def submit_job(
    bk_path: str, stripping_line: str, n_tuple_path: str, app, files_per_job=5
) -> None:
    """
    Submit a job to the grid, config defined in ./nTupleOptions.py

    The stripping line, bookkeeping path and nTuple path are passed across to this script

    Must provide an instantiated davinci application or something via app

    """
    check_bkfile_exists(bk_path)

    # Init job + tell Ganga where my config file is
    j = Job(name="Prompt DK3Pi Data")
    j.application = app
    j.application.options = ["./nTupleOptions.py"]
    j.application.platform = "x86_64-centos7-gcc8-opt"

    # Tell Ganga what data to use
    j.application.extraOpts = "print([i for i in range(10)])"

    # Job uses Dirac backend i guess (?)
    j.backend = Dirac()

    # Tell the job to split itself up a bit
    #j.splitter = SplitByFiles(filesPerJob=files_per_job)

    # Configure which files to spit out, i think?
    j.outputfiles = [DiracFile("*.root"), LocalFile("tmp.txt")]

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

    # Init DaVinci
    daVinci_app = create_davinci_application(".", "v45r1")

    # Create a job for one bk path to test
    year = "2011"
    submit_job(
        bookkeeping_paths[year][0],
        stripping_lines[year],
        "test.root",
        daVinci_app,
        1,
    )


# Unfortunately we can't wrap this in if name==main since we need to run it via ganga
main()
