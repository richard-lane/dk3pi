import os

# If the directory where I will be storing my local copy of DaVinci already exists, raise
if os.path.exists("./DaVinciDev_v45r1/"):
    raise Exception("rm the davnci dir before running")

# Create a local copy of DaVinci v45r1
j = Job(name='LHCb DK3Pi MC')
myApp = prepareGaudiExec('DaVinci','v45r1', myPath='.')

# Tell Ganga where my config file is
j.application = myApp
j.application.options = ['nTupleOptions.py']
j.application.platform = 'x86_64-centos7-gcc8-opt'

# Choose which MC to use
# Alternatives are specified in MCPaths.txt
bkPath = "/MC/2016/Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8/Sim09c/Trig0x6138160F/Reco16/Turbo03/Stripping28r1NoPrescalingFlagged/27265000/ALLSTREAMS.DST"
data  = BKQuery(bkPath, dqflag=['OK']).getDataset()

# Submit the job
j.inputdata = data
j.backend = Dirac()
j.splitter = SplitByFiles(filesPerJob=50)  # this seems like a good number to use
j.outputfiles = [DiracFile('*.root')]
j.submit()

