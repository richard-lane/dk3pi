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


def davinci_config(
    n_tuple_name: str,
    n_tuple_path: str,
    stripping_line: str,
    stream: str,
    decay_descriptor: str,
) -> str:
    """
    Return a str containing complete DaVinci config

    This is a kind of stupid way of doing things but it lets us manage DaVinci without having a separate script

        n_tuple_name: tree name or something within the nTuple. Basically just a unique identifier
        n_tuple_path: path to save our nTuple to
        stripping_line: stripping line from the STRIPPING project website
        stream: stream name; can also get this from the STRIPPING site. Probably "Charm" for real data
        decay_descriptor: decay descriptor suitable for passing to decayTreeTuple.setDescriptorTemplate. Must define a branch called "D" so we can find its decay time

    """
    # Return a string that will be written to a file that is used to configure DaVinci
    return (
        "print('Loading DaVinci config')\n"
        f"print('\tpath: {n_tuple_path}')\n"
        f"print('\tline: {stripping_line}')\n"
        f"print('\tstream: {stream}')\n"
        f"print('\tdecay: {decay_descriptor}')\n"
        "\n"
        "from Configurables import DecayTreeTuple\n"
        "from DecayTreeTuple.Configuration import *\n"
        "from Configurables import DaVinci\n"
        "from PhysConf.Filters import LoKi_Filters\n"
        "\n"
        "\n"
        "# We only want to look at a single TES location, so we don't need to unpack every event from the file\n"
        "# This means we should use a pre-filter\n"
        f"""DaVinci().EventPreFilters = LoKi_Filters(STRIP_Code = "HLT_PASS_RE('.*{stripping_line}.*')").filters('Filters')\n"""
        "\n"
        "# Initialise an nTuple for all of our data\n"
        f"dtt = DecayTreeTuple('{n_tuple_name}')\n"
        "\n"
        "# Find this location from bender dst-dump, or something\n"
        "# Since the data file is a microDST, this is relative to /Event/{stream}\n"
        f"dtt.Inputs = ['Phys/{stripping_line}/Particles']\n"
        "\n"
        f"dtt.setDescriptorTemplate('{decay_descriptor}')\n"
        "dtt.D.addTupleTool('TupleToolPropertime')\n"
        "\n"
        "# Configure DaVinci itself\n"
        "DaVinci().UserAlgorithms += [dtt]\n"
        "DaVinci().InputType = 'MDST'\n"
        f"DaVinci().TupleFile = '{n_tuple_path}'\n"
        "DaVinci().PrintFreq = 1000\n"
        "DaVinci().DataType = '2011'\n"
        "DaVinci().Simulation = False\n"
        "\n"
        "# This is necessary since we're reading from a microDST\n"
        f"DaVinci().RootInTES = '/Event/{stream}'\n"
        "\n"
        "# Ask for luminosity information\n"
        "DaVinci().Lumi = not DaVinci().Simulation\n"
        "DaVinci().CondDBtag = 'default'\n"
        "DaVinci().DDDBtag = 'default'\n"
        "\n"
        "# Ask for all events\n"
        "DaVinci().EvtMax = -1\n"
        "\n"
        "print('Loaded DaVinci config')\n"
    )


def submit_job(
    bk_path: str,
    stripping_line: str,
    n_tuple_path: str,
    n_tuple_name: str,
    stream: str,
    decay_descriptor: str,
    app,
    files_per_job=5,
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
    j.application.extraOpts = davinci_config(
        n_tuple_name, n_tuple_path, stripping_line, stream, decay_descriptor
    )

    # Job uses Dirac backend i guess (?)
    j.backend = Dirac()

    # Tell the job which data to use, and how to split the files up for processing
    j.inputdata = BKQuery(bk_path, dqflag=["OK"]).getDataset()
    j.splitter = SplitByFiles(filesPerJob=files_per_job)

    # Configure which files to spit out, i think?
    j.outputfiles = [DiracFile("*.root"), LocalFile("*.root")]

    # Submit job
    j.submit()


def main():
    # We want data from all years
    years = ("2011", "2012", "2015", "2016", "2017", "2018")
    bookkeeping_paths = bk_paths(years, ("MagDown", "MagUp"), "CHARM.MDST")

    # One stripping line per year that contains the data I want
    # From the LHCb stripping project site
    wrong_sign_prompt_descriptor = "${Dstar}[D*(2010)+ -> ${D}(D0 -> ${K}K+ ${pi1}pi- ${pi2}pi- ${pi3}pi+) ${pisoft}pi+]CC"
    right_sign_prompt_descriptor = "${Dstar}[D*(2010)+ -> ${D}(D0 -> ${K}K- ${pi1}pi+ ${pi2}pi+ ${pi3}pi-) ${pisoft}pi+]CC"
    prompt_stripping_lines = {
        "2011": "DstarPromptWithD02HHHHLine",
        "2012": "DstarPromptWithD02HHHHLine",
        "2015": "DstarD2HHHHDstarD2KPiPiPiLine ",
        "2016": "DstarD2HHHHDstarD2KPiPiPiLine ",
        "2017": "DstarD2HHHHDstarD2KPiPiPiLine",
        "2018": "DstarD2HHHHDstarD2KPiPiPiLine",
    }

    wrong_sign_semileptonic_descriptor = (
        "[B+ -> ${D}(D~0 -> ${K}K- ${pi1}pi+ ${pi2}pi+ ${pi3}pi-) ${mu}mu+]CC"
    )
    right_sign_semileptonic_descriptor = (
        "[B+ -> ${D}(D~0 -> ${K}K+ ${pi1}pi- ${pi2}pi- ${pi3}pi+) ${mu}mu+]CC"
    )
    semileptonic_stripping_lines = {
        "2011": "b2D0MuXK3PiCharmFromBSemiLine",
        "2012": "b2D0MuXK3PiCharmFromBSemiLine",
        "2015": "b2D0MuXK3PiCharmFromBSemiLine",
        "2016": "b2D0MuXK3PiCharmFromBSemiLine",
        "2017": "b2D0MuXK3PiCharmFromBSemiLine",
        "2018": "b2D0MuXK3PiCharmFromBSemiLine",
    }

    # Init DaVinci
    daVinci_app = create_davinci_application("tmp", "v45r1")

    # Create a job for one bk path to test
    submit_job(
        "/LHCb/Collision11/Beam3500GeV-VeloClosed-MagUp/Real Data/Reco14/Stripping21r1/90000000/CHARM.MDST",
        prompt_stripping_lines["2011"],
        "test_2011_WS_prompt.root",
        "dk3pi",
        "Charm",
        wrong_sign_prompt_descriptor,
        daVinci_app,
        5,
    )


# Unfortunately we can't wrap this in if name==main since we need to run it via ganga
main()
