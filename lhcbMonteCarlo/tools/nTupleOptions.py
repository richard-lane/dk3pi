from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import DaVinci

# Stream and stripping line
stream = "AllStreams"
line = "DstarD2HHHHDstarD2KPiPiPiLine"

# Create an nTuple
# Since we're dealing with an MDST the input location is relative to /Event/${stream}
dtt = DecayTreeTuple("TupleDstToD0pi_D0ToKpipipi")
dtt.Inputs = ["/Phys/{0}/Particles".format(line)]

# Add branches for each particle and the proper time of the D0
dtt.setDescriptorTemplate(
    "${Dstar}[D*(2010)+ -> ${D}(D0 -> ${K}K- ${pi1}pi+ ${pi2}pi+ ${pi3}pi-) ${pisoft}pi+]CC"
)
dtt.D.addTupleTool("TupleToolPropertime")

# Add some tuple tools that might be what I want
tuple_tools = {
    "TupleToolTISTOS",
    "TupleToolL0Calo",
    # "TupleToolTagging",
    # "TupleToolRecoStats",
    "TupleToolTrigger",
    "TupleToolMCTruth",  # MC
    "TupleToolMCBackgroundInfo",
    "TupleToolDalitz",  # MC
    # "TupleToolTrackInfo",
}
for tool in tuple_tools:
    dtt.addTupleTool(tool)

# Configure DaVinci
DaVinci().UserAlgorithms += [dtt]
DaVinci().InputType = "MDST"
DaVinci().TupleFile = "DVntuple.root"
DaVinci().PrintFreq = 1000
DaVinci().DataType = "2018"
DaVinci().Simulation = True

# Since we're reading from a microDST
DaVinci().RootInTES = "/Event/{0}".format(stream)

# Only ask for luminosity information when not using simulated data
DaVinci().Lumi = not DaVinci().Simulation
DaVinci().EvtMax = -1
DaVinci().CondDBtag = "sim-20190430-vc-md100"
DaVinci().DDDBtag = "dddb-20170721-3"

