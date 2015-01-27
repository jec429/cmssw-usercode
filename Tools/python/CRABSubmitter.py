#!/usr/bin/env python

import os,math,subprocess
import JChaves.Tools.Samples as Samples

def crab2_submit(*arg):
    # arguments are:
    # 0) name
    # 1) number of events
    # 2-) name modifiers
    crab_cfg = '''
[CRAB]
scheduler		= remoteGlidein
jobtype                 = cmssw
%(extra_CRAB)s

[CMSSW]
total_number_of_events  = %(nevents)i
events_per_job  = %(events_per)i
pset                    = %(pset)s
datasetpath             = %(data_set)s
output_file             = %(outfile)s
%(extra_CMSSW)s
    
[USER]
return_data             = 0
copy_data               = 1
storage_element         = T3_US_FNALLPC
ui_working_dir          = %(ui_dir)s
user_remote_dir         = %(sample)s
%(extra_USER)s  
    '''

    name_mod = ''
    if 'nstep' in arg[0]:
        if len(arg) == 3:
            name_mod = '_ee_'+str(arg[2])
        if len(arg) > 3:
            name_mod = '_'+arg[2]+'_'+str(arg[3])
    elif len(arg) > 2:
        name_mod = arg[2]
        
    datasets = ['None',Samples.ttbar_samples[0].dataset,'','']
    uidirs   = ['Nstep'+name_mod,'WR_TTbar_'+name_mod,'WR_DYJets_'+name_mod,'']
    samples  = ['Nstep'+name_mod,'WR_ttbar_histos_'+name_mod,'WR_dyjets_histos_'+name_mod,'']
    nevents =  arg[1]
    
    if 'nstep' in arg[0]:
        index = 0    
    elif 'ttbar' in arg[0]:
        index = 1
        if arg[1] == -1:
            nevents = Samples.ttbar_samples[0].size
    elif 'dyjets' in arg[0]:        
        index = 2
        assert len(arg) > 2
        if '100to200' in arg[2]:
            datasets[2] = Samples.dyjets_samples[0].dataset
            if arg[1] == -1:
                nevents = Samples.dyjets_samples[0].size
        elif '200to400' in arg[2]:
            datasets[2] = Samples.dyjets_samples[1].dataset
            if arg[1] == -1:
                nevents = Samples.dyjets_samples[1].size            
        elif '400to600' in arg[2]:
            datasets[2] = Samples.dyjets_samples[2].dataset
            if arg[1] == -1:
                nevents = Samples.dyjets_samples[2].size            
        elif '600toInf' in arg[2]:
            datasets[2] = Samples.dyjets_samples[3].dataset
            if arg[1] == -1:
                nevents = Samples.dyjets_samples[3].size            
            
    vd = locals()
    vd['extra_CRAB'] = ''
    vd['extra_CMSSW'] = ''
    vd['extra_USER'] = ''
    vd['nevents']     = nevents
    vd['events_per']  = math.ceil(nevents/500.0)
    vd['data_set']    = datasets[index]
    vd['ui_dir']      = 'crab/'+uidirs[index]
    vd['sample']      = samples[index]

    if index == 0:
        vd['pset'] = 'step1_cfg.py'
        vd['outfile'] = 'step1.root,step4_PAT.root'
        vd['extra_CMSSW'] = 'get_edm_output = 0\nignore_edm_output = 1'
        vd['extra_USER'] = 'additional_input_files  = step2_cfg.py,step3_cfg.py,step4_cfg.py\nscript_exe		= nstep.sh'
    elif 'genhistos' in name_mod:
        vd['pset']    = 'genhistos_cfg.py'
        vd['outfile'] = 'genhistos.root'        
    else:
        vd['pset']    = 'wr_analyzer_cfg.py'
        vd['outfile'] = 'out.root'
    
    open('crab.cfg','wt').write(crab_cfg % vd)
    os.system('crab -create')
    proc = subprocess.Popen(['tail','crab/'+uidirs[index]+'/share/arguments.xml'],stdout=subprocess.PIPE,)
    last_job = 499
    for line in proc.stdout:
        if 'JobID' in line:
            last_job = int(line.split(' ')[-2].split('"')[-2]) -1
            
    for i in range(0,last_job/500 + 1):
        print uidirs[index]
        os.system('crab -submit 500 -c crab/'+uidirs[index])
    os.system('mkdir -p psets')
    os.system('mv crab.cfg psets/crab_%s_cfg.py' % samples[index])

def crab3_submit(*arg):
    # arguments are:
    # 0) name
    # 1) number of events
    # 2-) name modifiers
    crab_cfg = '''
from CRABClient.client_utilities import getBasicConfig 
config = getBasicConfig()

config.General.requestName = '%(ui_dir)s'
config.General.workArea = 'crab'

config.JobType.pluginName = '%(plugin)s'
config.JobType.psetName = '%(pset)s'

config.Data.splitting = '%(split)s'
config.Data.unitsPerJob = %(nevents)i
NJOBS = %(njobs)i  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS

config.Data.outLFN = '/store/user/jchavesb/%(sample)s'
config.Site.storageSite = 'T3_US_FNALLPC'

'''
    name_mod = ''
    if 'nstep' in arg[0]:
        if len(arg) == 3:
            name_mod = '_ee_'+str(arg[2])
        if len(arg) > 3:
            name_mod = '_'+arg[2]+'_'+str(arg[3])
    elif len(arg) > 2:
        name_mod = arg[2]
        
    datasets = ['None',Samples.ttbar_samples[0].dataset,'','']
    uidirs   = ['Nstep'+name_mod,'WR_TTbar_'+name_mod,'WR_DYJets_'+name_mod,'WR_Signal_'+name_mod]
    samples  = ['Nstep'+name_mod,'WR_ttbar_histos_'+name_mod,'WR_dyjets_histos_'+name_mod,'WR_signal_histos_'+name_mod]
    nevents =  arg[1]
    
    if 'nstep' in arg[0]:
        index = 0
    elif 'ttbar' in arg[0]:
        index = 1
        if arg[1] == -1:
            #nevents = Samples.ttbar_samples[0].size
            nevents = Samples.ttbar_samples[0].nfiles
    elif 'dyjets' in arg[0]:        
        index = 2
        assert len(arg) > 2
        if '100to200' in arg[2]:
            datasets[2] = Samples.dyjets_samples[0].dataset
            if arg[1] == -1:
                nevents = Samples.ttbar_samples[0].nfiles
        elif '200to400' in arg[2]:
            datasets[2] = Samples.dyjets_samples[1].dataset
            if arg[1] == -1:
                nevents = Samples.ttbar_samples[0].nfiles
        elif '400to600' in arg[2]:
            datasets[2] = Samples.dyjets_samples[2].dataset
            if arg[1] == -1:
                nevents = Samples.ttbar_samples[0].nfiles
        elif '600toInf' in arg[2]:
            datasets[2] = Samples.dyjets_samples[3].dataset
            if arg[1] == -1:
                nevents = Samples.ttbar_samples[0].nfiles
    elif 'signal' in arg[0]:
        index = 3
        s = arg[0].split('_')[-1] #TODO: find the signal sample        
        datasets[3] = Samples.signal_gensamples[0].dataset
        nevents = Samples.signal_gensamples[0].nfiles
           
    vd = locals()
    vd['ui_dir'] = uidirs[index]
    vd['pset'] = ''
    vd['data_set'] = datasets[index]
    vd['nevents'] = math.ceil(nevents/500.0)
    vd['njobs'] = 500
    vd['sample']      = samples[index]

    if index == 0:
        vd['pset'] = 'step1_cfg.py'
        vd['split'] = 'EventBased'
        crab_cfg = crab_cfg + "config.JobType.scriptExe = 'nstep.sh'\n        config.JobType.inputFiles = ['step2_cfg.py','step3_cfg.py','step4_cfg.py']\n        config.JobType.outputFiles = ['step4.root']\n"
        #vd['outfile'] = 'step1.root,step4_PAT.root'
        #vd['extra_CMSSW'] = 'get_edm_output = 0\nignore_edm_output = 1'
        #vd['extra_USER'] = 'additional_input_files  = step2_cfg.py,step3_cfg.py,step4_cfg.py\nscript_exe		= nstep.sh'
    elif 'genhistos' in name_mod:
        vd['pset']    = 'genhistos_cfg.py'
        vd['plugin']  = 'Analysis'
        vd['split'] = 'FileBased'
        crab_cfg = crab_cfg + "config.Data.inputDataset = '" + datasets[index] + "'\n"
    else:
        vd['pset']    = 'wr_analyzer_cfg.py'
        vd['plugin']  = 'Analysis'
        vd['split'] = 'FileBased'
        crab_cfg = crab_cfg + "config.Data.inputDataset = '" + datasets[index] + "'\n"
    
    open('crabConfig.py','wt').write(crab_cfg % vd)
    os.system('crab submit')
    os.system('mkdir -p psets')
    os.system('mv crabConfig.py psets/crabConfig_%s.py' % samples[index])
    
