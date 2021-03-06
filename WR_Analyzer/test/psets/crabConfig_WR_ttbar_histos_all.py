
from CRABClient.client_utilities import getBasicConfig 
config = getBasicConfig()

config.General.requestName = 'WR_TTbar_all'
config.General.workArea = 'crab'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'wr_analyzer_cfg.py'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
NJOBS = 500  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS

config.Data.outLFN = '/store/user/jchavesb/WR_ttbar_histos_all'
config.Site.storageSite = 'T3_US_FNALLPC'

config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
