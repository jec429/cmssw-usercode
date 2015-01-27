
from CRABClient.client_utilities import getBasicConfig 
config = getBasicConfig()

config.General.requestName = 'WR_Signal_genhistos'
config.General.workArea = 'crab'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'genhistos_cfg.py'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 100
NJOBS = 500  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS

config.Data.outLFN = '/store/user/jchavesb/WR_signal_histos_genhistos'
config.Site.storageSite = 'T3_US_FNALLPC'

config.Data.inputDataset = '/WRToNuEToEEJJ_MW-2000_MNu-1000_TuneCUETP8M1_13TeV-pythia8/RunIIFall14GS-MCRUN2_71_V1-v1/GEN-SIM'
