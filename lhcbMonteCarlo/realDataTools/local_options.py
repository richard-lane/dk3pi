"""
Local test the DaVinci config on an example microDST

Find the dtt input and stream from bender dst-dump

"""
from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import DaVinci
from GaudiConf import IOHelper
from PhysConf.Filters import LoKi_Filters

# We only want to look at a single TES location, so we don't need to unpack every event from the file
# This means we can use a pre-filter
DaVinci().EventPreFilters = LoKi_Filters(STRIP_Code = "HLT_PASS_RE('.*b2D0MuXK3PiCharmFromBSemiLine.*')").filters('Filters')

# Stripping line from stripping project website; remove the word 'Stripping' from the front
line = 'b2D0MuXK3PiCharmFromBSemiLine'
stream = 'Charm'

# Initialise an nTuple for all of our data
dtt = DecayTreeTuple('TupleDstToD0pi_D0ToKpipipi')

# Find this location from bender dst-dump, or something
# Since the data file is a microDST, this is relative to /Event/{stream}
dtt.Inputs = ['Phys/{0}/Particles'.format(line)]
dtt.setDescriptorTemplate('[B+ -> ${D}(D~0 -> ${K}K+ ${pi1}pi- ${pi2}pi- ${pi3}pi+) ${mu}mu+]CC')

# Add the proper decay time of the D0\n"
dtt.D.addTupleTool('TupleToolPropertime')

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
    './00041838_00000067_1.charm.mdst'
], clear=True)

