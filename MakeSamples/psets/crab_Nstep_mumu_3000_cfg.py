
[CRAB]
scheduler		= remoteGlidein
jobtype                 = cmssw


[CMSSW]
total_number_of_events  = 100000
events_per_job  = 200
pset                    = step1_cfg.py
datasetpath             = None
output_file             = step1.root,step4_PAT.root
get_edm_output = 0
ignore_edm_output = 1
    
[USER]
return_data             = 0
copy_data               = 1
storage_element         = T3_US_FNALLPC
ui_working_dir          = crab/Nstep_mumu_3000
user_remote_dir         = Nstep_mumu_3000
additional_input_files  = step2_cfg.py,step3_cfg.py,step4_cfg.py
script_exe		= nstep.sh  
    