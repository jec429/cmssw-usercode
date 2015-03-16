#!/usr/bin/env python

class Sample:
    def __init__(self,name,xsection,dataset,size,nfiles):
        self.xsection = xsection
        self.name = name
        self.dataset = dataset
        self.size = size
        self.nfiles = nfiles

################################################################################
################################################################################
# MINIAOD Samples:
################################################################################
################################################################################

signal_samples = [Sample('WR_ee_800',3.65E+00,'/eos/uscms/store/user/jchaves/Nstep_ee/',1E5,1),Sample('WR_ee_1000',1.78E+00,'/eos/uscms/store/user/jchaves/Nstep_ee_1000/',1E5,1),
                  Sample('WR_ee_1200',6.63E-01,'/eos/uscms/store/user/jchaves/Nstep_ee_1200/',1E5,1),Sample('WR_ee_1400',3.89E-01,'/eos/uscms/store/user/jchaves/Nstep_ee_1400/',1E5,1),
                  Sample('WR_ee_1600',1.77E-01,'/eos/uscms/store/user/jchaves/Nstep_ee_1600/',1E5,1),Sample('WR_ee_1800',1.17E-01,'/eos/uscms/store/user/jchaves/Nstep_ee_1800/',1E5,1),
                  Sample('WR_ee_2000',7.07E-02,'/eos/uscms/store/user/jchaves/Nstep_ee_2000/',1E5,1),Sample('WR_ee_2200',4.50E-02,'/eos/uscms/store/user/jchaves/Nstep_ee_2200/',1E5,1),
                  Sample('WR_ee_2400',2.48E-02,'/eos/uscms/store/user/jchaves/Nstep_ee_2400/',1E5,1),Sample('WR_ee_2600',1.50E-02,'/eos/uscms/store/user/jchaves/Nstep_ee_2600/',1E5,1),
                  Sample('WR_ee_2800',9.13E-03,'/eos/uscms/store/user/jchaves/Nstep_ee_2800/',1E5,1),Sample('WR_ee_3000',5.76E-03,'/eos/uscms/store/user/jchaves/Nstep_ee_3000/',1E5,1),
                  Sample('WR_ee_3200',3.40E-03,'/eos/uscms/store/user/jchaves/Nstep_ee_3200/',1E5,1),Sample('WR_ee_3400',2.63E-03,'/eos/uscms/store/user/jchaves/Nstep_ee_3400/',1E5,1),
                  Sample('WR_ee_3600',1.54E-03,'/eos/uscms/store/user/jchaves/Nstep_ee_3600/',1E5,1),Sample('WR_ee_3800',1.19E-03,'/eos/uscms/store/user/jchaves/Nstep_ee_3800/',1E5,1),
                  Sample('WR_ee_4000',8.01E-04,'/eos/uscms/store/user/jchaves/Nstep_ee_4000/',1E5,1),Sample('WR_ee_4200',4.73E-04,'/eos/uscms/store/user/jchaves/Nstep_ee_4200/',1E5,1),
                  Sample('WR_ee_4400',3.75E-04,'/eos/uscms/store/user/jchaves/Nstep_ee_4400/',1E5,1),Sample('WR_ee_4600',1.90E-04,'/eos/uscms/store/user/jchaves/Nstep_ee_4600/',1E5,1),
                  Sample('WR_ee_4800',1.52E-04,'/eos/uscms/store/user/jchaves/Nstep_ee_4800/',1E5,1),Sample('WR_ee_5000',9.12E-05,'/eos/uscms/store/user/jchaves/Nstep_ee_5000/',1E5,1),
                  Sample('WR_ee_5200',6.65E-05,'/eos/uscms/store/user/jchaves/Nstep_ee_5200/',1E5,1),Sample('WR_ee_5400',4.49E-05,'/eos/uscms/store/user/jchaves/Nstep_ee_5400/',1E5,1),
                  Sample('WR_ee_5600',2.54E-05,'/eos/uscms/store/user/jchaves/Nstep_ee_5600/',1E5,1),Sample('WR_ee_5800',2.02E-05,'/eos/uscms/store/user/jchaves/Nstep_ee_5800/',1E5,1),
                  Sample('WR_ee_6000',1.44E-05,'/eos/uscms/store/user/jchaves/Nstep_ee_6000/',1E5,1)]        

ttbar_samples = [Sample('ttbar',816,'/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',25446993,552)]

dyjets_samples = [Sample('dyjets_HT100to200',194.3,'/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',4054159,98),
                  Sample('dyjets_HT200to400',50.24,'/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',4666496,105),
                  Sample('dyjets_HT400to600',6.546,'/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',4931372,115),
                  Sample('dyjets_HT600toInf',2.179,'/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',4493574,108),
                  
                  Sample('dyjets_M200to400',2008.4,'/DYJetsToEEMuMu_M-200To400_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',26430,1),
                  Sample('dyjets_M400to800',2008.4,'/DYJetsToEEMuMu_M-400To800_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',48205,1),
                  Sample('dyjets_M800to1400',2008.4,'/DYJetsToEEMuMu_M-800To1400_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',24003,1),
                  Sample('dyjets_M1400to2300',2008.4,'/DYJetsToEEMuMu_M-1400To2300_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',23204,1),
                  Sample('dyjets_M3500to4500',2008.4,'/DYJetsToEEMuMu_M-3500To4500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',23297,1),
                  Sample('dyjets_M4500to6000',2008.4,'/DYJetsToEEMuMu_M-4500To6000_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',25925,1),
                  Sample('dyjets_M6000to7500',2008.4,'/DYJetsToEEMuMu_M-6000To7500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',28710,1),
                  Sample('dyjets_M7500to8500',2008.4,'/DYJetsToEEMuMu_M-7500To8500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',30151,2)
                  ]

EE_signal_gensamples = [Sample('WR_EE_2000',7.07E-02,'/WRToNuEToEEJJ_MW-2000_MNu-1000_TuneCUETP8M1_13TeV-pythia8/RunIIFall14GS-MCRUN2_71_V1-v1/GEN-SIM',50000,32),
                     Sample('WR_EE_5400',4.49E-05,'/WRToNuEToEEJJ_MW-5400_MNu-2700_TuneCUETP8M1_13TeV-pythia8/RunIIFall14GS-MCRUN2_71_V1-v1/GEN-SIM',49764,1)]
    
MUMU_signal_gensamples = [Sample('WR_MUMU_2000',6.22E-02,'/WRToNuMuToMuMuJJ_MW-2000_MNu-1000_TuneCUETP8M1_13TeV-pythia8/RunIIFall14GS-MCRUN2_71_V1-v1/GEN-SIM',50000,25)]
