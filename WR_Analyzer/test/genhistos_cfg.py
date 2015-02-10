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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring('file:/uscms/home/jchaves/private/CMSSW_7_2_2_patch2/src/MakeSamples/step0.root'),
                            )
# This only works interactively
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
                        
#process.source.fileNames = file_list(dset+'step1_*.root',True)
#process.source.fileNames = file_list('/WRToNuEToEEJJ_MW-2000_MNu-1000_TuneCUETP8M1_13TeV-pythia8/RunIIFall14GS-MCRUN2_71_V1-v1/GEN-SIM',False)
outfile = 'genhistos0.root'

if ttbar:
    process.source.fileNames = file_list(Samples.ttbar_samples[0].dataset,False)
    outfile = 'genhistos_ttbar.root'
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

if dyjets:
    process.source.fileNames = file_list(Samples.dyjets_samples[0].dataset,False)
    outfile = 'genhistos_dyjets.root'
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.TFileService = cms.Service('TFileService', fileName = cms.string(outfile))


process.ana = cms.EDAnalyzer('GenHistos',
                            gen_src = cms.InputTag('genParticles'),
                            gen_jet_src = cms.InputTag('ak4GenJets'),
                            jet_pt_min = cms.double(40),
                            jet_eta_max = cms.double(2.5),
                            dilepton_mass_cut = cms.double(0.0),
                            lljj_mass_cut = cms.double(0.0)
                            )
process.ana2 = process.ana.clone(dilepton_mass_cut = cms.double(200.0))
process.ana3 = process.ana.clone(dilepton_mass_cut = cms.double(200.0),
                                 lljj_mass_cut = cms.double(600.0))

process.p = cms.Path(process.ana)
process.p2 = cms.Path(process.ana2)
process.p3 = cms.Path(process.ana3)

###################################################################################################
###################################################################################################
# Use my CRABSubmitter for now:
# crab2_submit(name,nevents,name modifiers)
# name = nstep -> MakeSamples
# name = ttbar
# name = dyjets
if __name__ == '__main__' and hasattr(sys, 'argv') and submit:
    from JChaves.Tools.CRABSubmitter import *
    crab3_submit('signal_2000',-1,'genhistos')
    if ttbar:
        #crab2_submit('ttbar',-1,'genhistos')
        crab3_submit('ttbar',-1,'genhistos')
    if dyjets:
        for x in ['HT-100to200_genhistos','HT-200to400_genhistos','HT-400to600_genhistos','HT-600toInf_genhistos']:
            crab2_submit('dyjets',-1,x)
