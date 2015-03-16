from JChaves.Tools.CMSTools import *
import JChaves.Tools.Samples as Samples
import os,glob,sys
import FWCore.ParameterSet.Config as cms

ttbar = 'ttbar' in sys.argv
dyjets = 'dyjets' in sys.argv
submit = 'submit' in sys.argv

process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring('root://cmsxrootd-site.fnal.gov//store/mc/Summer12/WRToNuLeptonToLLJJ_MW-3500_MNu-1750_TuneZ2star_8TeV-pythia6-tauola/GEN-SIM/START50_V13-v2/0000/04317A88-1FA2-E111-8D7A-485B39897242.root'),
                            )
# This only works interactively
mass = '_'
if not ttbar and not dyjets:
    for m in sys.argv:
        if '00' in m:
            mass = '_'+m
dset = ''
for s in Samples.signal_samples:
    if mass in s.name:
        dset = s.dataset
        
#process.source.fileNames = file_list(dset+'step4_PAT*.root',True)
process.source.fileNames = file_list('/eos/uscms/store/user/jchaves/Nstep_MUMU_2000_reco/EXO-Phys14DR-00009_*.root',True)
#process.source.fileNames = file_list('/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',False)
outfile = 'xsec.root'

if ttbar:
    process.source.fileNames = file_list(Samples.ttbar_samples[0].dataset,False)
    outfile = 'xsec_ttbar.root'
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

if dyjets:
    process.source.fileNames = file_list(Samples.dyjets_samples[0].dataset,False)
    outfile = 'xsec_dyjets.root'
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.TFileService = cms.Service('TFileService', fileName = cms.string(outfile))

process.load('JChaves.WR_Analyzer.xsec_cfi')


process.ana.scale = cms.double(194.3)
#process.ana.scale = cms.double(50.24)
#process.ana.scale = cms.double(6.546)
#process.ana.scale = cms.double(2.179)

process.p = cms.Path(process.ana)


###################################################################################################
###################################################################################################
# Use my CRABSubmitter for now:
# crab2_submit(name,nevents,name modifiers)
# name = nstep -> MakeSamples
# name = ttbar
# name = dyjets
if __name__ == '__main__' and hasattr(sys, 'argv') and submit:
    from JChaves.Tools.CRABSubmitter import *
    
    if ttbar:
        #crab2_submit('ttbar',-1,'all')
        crab3_submit('ttbar',-1,'all')
    if dyjets:
        for x in ['HT-100To200','HT-200To400','HT-400To600','HT-600ToInf']:
        #for x in ['M-200To400','M-400To800','M-800To1400','M-1400To2300','M-3500To4500','M-4500To6000','M-6000To7500','M-7500To8500',]:
            #crab2_submit('dyjets',-1,x)
            crab3_submit('dyjets',-1,x)
