
[CRAB]
scheduler		= remoteGlidein
jobtype                 = cmssw


[CMSSW]
total_number_of_events  = 25446993
events_per_job  = 50894
pset                    = genhistos_cfg.py
datasetpath             = /TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
output_file             = genhistos.root

    
[USER]
return_data             = 0
copy_data               = 1
storage_element         = T3_US_FNALLPC
ui_working_dir          = crab/WR_TTbar_genhistos
user_remote_dir         = WR_ttbar_histos_genhistos
  
    