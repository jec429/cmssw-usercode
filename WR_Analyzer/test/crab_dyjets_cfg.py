
[CMSSW]
total_number_of_events  = 5000000
events_per_job  = 10021
pset                    = wr_analyzer_cfg.py
datasetpath             = /DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
output_file             = out.root
    
[USER]
return_data             = 0
copy_data               = 1
storage_element         = T3_US_FNALLPC
ui_working_dir          = crab/WR_DYJets
user_remote_dir         = WR_dyjets_histos
    
    
[CRAB]
scheduler		= remoteGlidein
jobtype                 = cmssw
    