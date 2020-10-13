"""
Local test the DaVinci config on an example microDST

Find the dtt input and stream from bender dst-dump

"""
from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import DaVinci
from GaudiConf import IOHelper

# Stripping line from stripping project website; remove the word 'Stripping' from the front
line = 'DstarPromptWithD02HHHHLine'
stream = 'Charm'

# Initialise an nTuple for all of our data
dtt = DecayTreeTuple('TupleDstToD0pi_D0ToKpipipi')

# Find this location from bender dst-dump, or something
# Since the data file is a microDST, this is relative to /Event/{stream}
dtt.Inputs = ['Phys/{0}/Particles'.format(line)]
dtt.Decay = '[D*(2010)+ -> ^(D0 -> ^K- ^pi+ ^pi+ ^pi-) ^pi+]CC'

# Add some tuple tools\n"
tuple_tools = {
    'TupleToolTISTOS',
    'TupleToolL0Calo',
    'TupleToolRecoStats',
    'TupleToolTrigger',
    'TupleToolTrackInfo',
}
for tool in tuple_tools:
    dtt.addTupleTool(tool)

# Add branches for each particle that we're interested in
dtt.addBranches(
    {
        'Dstar': '[D*(2010)+ -> (D0 -> K- pi+ pi+ pi-) pi+]CC',
        'D0': '[D*(2010)+ -> ^(D0 -> K- pi+ pi+ pi-) pi+]CC',
        'Kminus': '[D*(2010)+ -> (D0 -> ^K- pi+ pi+ pi-) pi+]CC',
        'pi1plus': '[D*(2010)+ -> (D0 -> K- ^pi+ pi+ pi-) pi+]CC',
        'pi2plus': '[D*(2010)+ -> (D0 -> K- pi+ ^pi+ pi-) pi+]CC',
        'pi3minus': '[D*(2010)+ -> (D0 -> K- pi+ pi+ ^pi-) pi+]CC',
        'pisoft': '[D*(2010)+ -> (D0 -> K- pi+ pi+ pi-) ^pi+]CC',
    }
)

# Add the proper decay time of the D0\n"
dtt.D0.addTupleTool('TupleToolPropertime')

# Configure DaVinci itself
DaVinci().UserAlgorithms += [dtt]
DaVinci().InputType = 'MDST'
DaVinci().TupleFile = 'test_local.root'
DaVinci().PrintFreq = 1000
DaVinci().DataType = '2011'
DaVinci().Simulation = False

# This is necessary since we're reading from a microDST
DaVinci().RootInTES = '/Event/{0}'.format(stream)

# Ask for luminosity information
DaVinci().Lumi = not DaVinci().Simulation
DaVinci().CondDBtag = 'default'
DaVinci().DDDBtag = 'default'

# Ask for all events
DaVinci().EvtMax = -1

# Use local input data
IOHelper().inputFiles([
    './00041838_00000057_1.charm.mdst'
], clear=True)

