
[CRAB]
scheduler		= remoteGlidein
jobtype                 = cmssw


[CMSSW]
total_number_of_events  = 4493574
events_per_job  = 8988
pset                    = genhistos_cfg.py
datasetpath             = /DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
output_file             = genhistos.root

    
[USER]
return_data             = 0
copy_data               = 1
storage_element         = T3_US_FNALLPC
ui_working_dir          = crab/WR_DYJets_HT-600toInf_genhistos
user_remote_dir         = WR_dyjets_histos_HT-600toInf_genhistos
  
    