import sys

from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import DaVinci

# Arg parsing
# Arguments are passed through from ganga to this conf file via sys.argv, i think
# Should be a list of [bk path, stripping line, root file path]
print("Hello from the davinci script")
with open("tmp.txt", "w") as f:
    f.write("writing sys argv:\n")
    f.write(str(sys.argv))
    f.write("\nsys argv written\n")

## Initialise an nTuple for all of our data
#dtt = DecayTreeTuple("TupleDstToD0pi_D0ToKpipipi")
#dtt.Inputs = dst_paths
#dtt.Decay = "[D*(2010)+ -> ^(D0 -> ^K- ^pi+ ^pi+ ^pi-) ^pi+]CC"
#
## Add some tuple tools
#tuple_tools = {
#    "TupleToolTISTOS",
#    "TupleToolL0Calo",
#    "TupleToolRecoStats",
#    "TupleToolTrigger",
#    "TupleToolTrackInfo",
#}
#for tool in tuple_tools:
#    dtt.addTupleTool(tool)
#
## Add branches for each particle that we're interested in
#dtt.addBranches(
#    {
#        "Dstar": "[D*(2010)+ -> (D0 -> K- pi+ pi+ pi-) pi+]CC",
#        "D0": "[D*(2010)+ -> ^(D0 -> K- pi+ pi+ pi-) pi+]CC",
#        "Kminus": "[D*(2010)+ -> (D0 -> ^K- pi+ pi+ pi-) pi+]CC",
#        "pi1plus": "[D*(2010)+ -> (D0 -> K- ^pi+ pi+ pi-) pi+]CC",
#        "pi2plus": "[D*(2010)+ -> (D0 -> K- pi+ ^pi+ pi-) pi+]CC",
#        "pi3minus": "[D*(2010)+ -> (D0 -> K- pi+ pi+ ^pi-) pi+]CC",
#        "pisoft": "[D*(2010)+ -> (D0 -> K- pi+ pi+ pi-) ^pi+]CC",
#    }
#)
#
## Add the proper decay time of the D0
#dtt.D0.addTupleTool("TupleToolPropertime")
#
## Configure DaVinci itself
#DaVinci().UserAlgorithms += [dtt]
#DaVinci().InputType = "DST"
#DaVinci().TupleFile = "prompt.root"
#DaVinci().PrintFreq = 1000
#
## Ask for luminosity information
#DaVinci().Lumi = not DaVinci().Simulation
#
## Ask for all events
#DaVinci().EvtMax = -1
