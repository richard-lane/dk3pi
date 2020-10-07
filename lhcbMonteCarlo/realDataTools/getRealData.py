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


def create_davinci_application(path: str, version: str):
    """
    Create a local copy of the Davinic application

    returns whatever prepareGaudiExec() returns, idk i can't find any documentation on it ???

    """
    if os.path.exists(os.path.join(path, "DaVinciDev_v45r1/")):
        raise Exception("rm the davnci dir before submitting job")

    return prepareGaudiExec("DaVinci", version, myPath=path)


def davinci_config(n_tuple_path: str, stripping_line: str) -> str:
    """
    Return a str containing complete DaVinci config

    This is a kind of stupid way of doing things but it lets us manage DaVinci without having a separate script

    """
    return "print('Loading DaVinci config')"
    "from Configurables import DecayTreeTuple\n"
    "from DecayTreeTuple.Configuration import *\n"
    "from Configurables import DaVinci\n"
    "\n"
    "\n"
    "# Initialise an nTuple for all of our data\n"
    "dtt = DecayTreeTuple('TupleDstToD0pi_D0ToKpipipi')\n"
    f"dtt.Inputs = ['/Event/AllStreams/Phys/{stripping_line}/Particles']\n"
    "dtt.Decay = '[D*(2010)+ -> ^(D0 -> ^K- ^pi+ ^pi+ ^pi-) ^pi+]CC'\n"
    "\n"
    "# Add some tuple tools\n"
    "#tuple_tools = {\n"
    "#    'TupleToolTISTOS',\n"
    "#    'TupleToolL0Calo',\n"
    "#    'TupleToolRecoStats',\n"
    "#    'TupleToolTrigger',\n"
    "#    'TupleToolTrackInfo',\n"
    "#}\n"
    "#for tool in tuple_tools:\n"
    "#    dtt.addTupleTool(tool)\n"
    "\n"
    "# Add branches for each particle that we're interested in\n"
    "dtt.addBranches(\n"
    "    {\n"
    "        'Dstar': '[D*(2010)+ -> (D0 -> K- pi+ pi+ pi-) pi+]CC',\n"
    "        'D0': '[D*(2010)+ -> ^(D0 -> K- pi+ pi+ pi-) pi+]CC',\n"
    "        'Kminus': '[D*(2010)+ -> (D0 -> ^K- pi+ pi+ pi-) pi+]CC',\n"
    "        'pi1plus': '[D*(2010)+ -> (D0 -> K- ^pi+ pi+ pi-) pi+]CC',\n"
    "        'pi2plus': '[D*(2010)+ -> (D0 -> K- pi+ ^pi+ pi-) pi+]CC',\n"
    "        'pi3minus': '[D*(2010)+ -> (D0 -> K- pi+ pi+ ^pi-) pi+]CC',\n"
    "        'pisoft': '[D*(2010)+ -> (D0 -> K- pi+ pi+ pi-) ^pi+]CC',\n"
    "    }\n"
    ")\n"
    "\n"
    "# Add the proper decay time of the D0\n"
    "dtt.D0.addTupleTool('TupleToolPropertime')\n"
    "\n"
    "# Configure DaVinci itself\n"
    "DaVinci().UserAlgorithms += [dtt]\n"
    "DaVinci().InputType = 'DST'\n"
    f"DaVinci().TupleFile = '{n_tuple_path}'\n"
    "DaVinci().PrintFreq = 1000\n"
    "\n"
    "# Ask for luminosity information\n"
    "DaVinci().Lumi = not DaVinci().Simulation\n"
    "\n"
    "# Ask for all events\n"
    "DaVinci().EvtMax = -1\n"


def submit_job(
    bk_path: str, stripping_line: str, n_tuple_path: str, app, files_per_job=5
) -> None:
    """
    Submit a job to the grid, config defined in ./nTupleOptions.py

    The stripping line, bookkeeping path and nTuple path are passed across to this script

    Must provide an instantiated davinci application or something via app

    """
    # Init job
    j = Job(name=f"Prompt DK3Pi Data: {n_tuple_path}")
    j.application = app
    j.application.platform = "x86_64-centos7-gcc8-opt"

    # Provide options to DaVinci via this string
    j.application.extraOpts = davinci_config(n_tuple_path, stripping_line)

    # Job uses Dirac backend i guess (?)
    j.backend = Dirac()

    # Tell the job which data to use, and how to split the files up for processing
    j.inputdata = BKQuery(bk_path, dqflag=["OK"]).getDataset()
    j.splitter = SplitByFiles(filesPerJob=files_per_job)

    # Configure which files to spit out, i think?
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

    # Init DaVinci
    daVinci_app = create_davinci_application(".", "v45r1")

    # Create a job for one bk path to test
    year = "2011"
    submit_job(
        bookkeeping_paths[year][0], stripping_lines[year], "test.root", daVinci_app, 1
    )


# Unfortunately we can't wrap this in if name==main since we need to run it via ganga
main()
