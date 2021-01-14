import os

# If the directory where I will be storing my local copy of DaVinci already exists, raise
if os.path.exists("./DaVinciDev_v45r1/"):
    raise Exception("rm the davnci dir before running")

# Create a local copy of DaVinci v45r1
j = Job(name='2018 WS MC')
myApp = prepareGaudiExec('DaVinci','v45r1', myPath='.')

# Tell Ganga where my config file is
j.application = myApp
j.application.options = ['nTupleOptions.py']
j.application.platform = 'x86_64-centos7-gcc8-opt'

# Choose which MC to use
# Alternatives are specified in MCPaths.txt
bkPath = "/MC/2018/Beam6500GeV-2018-MagDown-Nu1.6-25ns-Pythia8/Sim09j-ReDecay01/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34NoPrescalingFlagged/Turbo05Filtered/27165072/DSTARD02HHHH.HLTFILTER.MDST"
data  = BKQuery(bkPath, dqflag=['OK']).getDataset()

# Submit the job
j.inputdata = data
j.backend = Dirac()
j.splitter = SplitByFiles(filesPerJob=1)
j.outputfiles = [LocalFile('*.root')]
j.submit()
