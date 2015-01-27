from CRABClient.client_utilities import getBasicConfig 
config = getBasicConfig()

config.General.requestName = 'Nstep_MuMu_10'
config.General.workArea = 'crab/Nstep_MuMu_10'

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'step1_cfg.py'
config.JobType.scriptExe = 'nstep.sh'
config.JobType.inputFiles = ['step2_cfg.py','step3_cfg.py','step4_cfg.py']
config.JobType.outputFiles = ['step4.root']

#config.Data.primaryDataset = 'MinBias'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 10
NJOBS = 1  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
#config.Data.outLFN = '/store/user/jchaves/Nstep_MuMu_1000/' # or '/store/group/<subdir>'
#config.Data.publication = True
#config.Data.publishDataName = 'CRAB3_tutorial_MC_generation_test2'

config.Site.storageSite = 'T3_US_FNALLPC'
