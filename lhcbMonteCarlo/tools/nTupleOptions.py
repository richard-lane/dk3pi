from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import DaVinci
from GaudiConf import IOHelper

# Stream and stripping line
stream = "AllStreams"
line = "DstarD2HHHHDstarD2KPiPiPiLine"

# Create an nTuple
dtt = DecayTreeTuple("TupleDstToD0pi_D0ToKpipipi")
dtt.Inputs = ["/Event/{0}/Phys/{1}/Particles".format(stream, line)]
dtt.Decay = "[D*(2010)+ -> ^(D0 -> ^K- ^pi+ ^pi+ ^pi-) ^pi+]CC"

# Add some tuple tools that might be what I want
tuple_tools ={'TupleToolTISTOS',
              'TupleToolL0Calo',
              # 'TupleToolTagging',
              'TupleToolRecoStats',
              'TupleToolTrigger',
              'TupleToolMCTruth', # MC
              'TupleToolMCBackgroundInfo',
              'TupleToolDalitz', # MC
              'TupleToolTrackInfo'}

for tool in tuple_tools:
    dtt.addTupleTool(tool)

# Add branches for each particle and the proper time of the D0
dtt.addBranches({"Dstar"    : "[D*(2010)+ -> (D0 -> K- pi+ pi+ pi-) pi+]CC",
                 "D0"       : "[D*(2010)+ -> ^(D0 -> K- pi+ pi+ pi-) pi+]CC",
                 "Kminus"   : "[D*(2010)+ -> (D0 -> ^K- pi+ pi+ pi-) pi+]CC",
                 "pi1plus"  : "[D*(2010)+ -> (D0 -> K- ^pi+ pi+ pi-) pi+]CC",
                 "pi2plus"  : "[D*(2010)+ -> (D0 -> K- pi+ ^pi+ pi-) pi+]CC",
                 "pi3minus" : "[D*(2010)+ -> (D0 -> K- pi+ pi+ ^pi-) pi+]CC",
                 "pisoft"   : "[D*(2010)+ -> (D0 -> K- pi+ pi+ pi-) ^pi+]CC"})
dtt.D0.addTupleTool("TupleToolPropertime")

# Configure DaVinci
DaVinci().UserAlgorithms += [dtt]
DaVinci().InputType = 'DST'
DaVinci().TupleFile = 'DVntuple.root'
DaVinci().PrintFreq = 1000
DaVinci().DataType = '2016'
DaVinci().Simulation = True

# Only ask for luminosity information when not using simulated data
DaVinci().Lumi = not DaVinci().Simulation
DaVinci().EvtMax = -1
DaVinci().CondDBtag = "sim-20170721-2-vc-md100"
DaVinci().DDDBtag = 'dddb-20170721-3'

