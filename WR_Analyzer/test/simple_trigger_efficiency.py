from JChaves.Tools.CMSTools import *
import JChaves.Tools.Samples as Samples
import os,glob,sys
import FWCore.ParameterSet.Config as cms

process = cms.Process('BasicAnalyzer')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(50000))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:pat.root'))
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source.fileNames = file_list('/eos/uscms/store/user/jchaves/Nstep_EE_2000_reco/EXO-Phys14DR-00009_*.root',True)
outfile = 'simple_trigger_efficiency.root'
process.TFileService = cms.Service('TFileService', fileName = cms.string(outfile))

process.genMus = cms.EDFilter('CandViewSelector', src = cms.InputTag('genParticles'), cut = cms.string('abs(pdgId) == 13 && abs(mother.pdgId) == 24'))
process.genMuCount = cms.EDFilter('CandViewCountFilter', src = cms.InputTag('genMus'), minNumber = cms.uint32(1))
                                
process.genEls = cms.EDFilter('CandViewSelector', src = cms.InputTag('genParticles'), cut = cms.string('abs(pdgId) == 11 && abs(mother.pdgId) == 24'))
process.genElCount = cms.EDFilter('CandViewCountFilter', src = cms.InputTag('genEls'), minNumber = cms.uint32(1))
                                
process.genMusInAcc = cms.EDFilter('CandViewSelector', src = cms.InputTag('genParticles'), cut = cms.string('abs(pdgId) == 13 && abs(mother.pdgId) == 24 && pt > 26 && abs(eta) < 2.1'))
process.genElsInAcc = cms.EDFilter('CandViewSelector', src = cms.InputTag('genParticles'), cut = cms.string('abs(pdgId) == 11 && abs(mother.pdgId) == 24 && pt > 30 && abs(eta) < 2.5'))
process.genMuInAccCount = cms.EDFilter('CandViewCountFilter', src = cms.InputTag('genMusInAcc'), minNumber = cms.uint32(1))
process.genElInAccCount = cms.EDFilter('CandViewCountFilter', src = cms.InputTag('genElsInAcc'), minNumber = cms.uint32(1))

process.RandomNumberGeneratorService = cms.Service('RandomNumberGeneratorService')
process.RandomNumberGeneratorService.SimpleTriggerEfficiency = cms.PSet(initialSeed = cms.untracked.uint32(1219))
process.RandomNumberGeneratorService.SimpleTriggerEfficiencyMu = cms.PSet(initialSeed = cms.untracked.uint32(1220))
process.RandomNumberGeneratorService.SimpleTriggerEfficiencyMuInAcc = cms.PSet(initialSeed = cms.untracked.uint32(1221))

#import prescales
process.SimpleTriggerEfficiency = cms.EDAnalyzer('SimpleTriggerEfficiency',
                                         trigger_results_src = cms.InputTag('TriggerResults', '', 'HLT'),
                                         weight_src = cms.InputTag(''),
                                         )

def setup_endpath(process, weight_src = ''):
    process.SimpleTriggerEfficiency = SimpleTriggerEfficiency.clone(weight_src = weight_src)
    process.SimpleTriggerEfficiency.trigger_results_src = cms.InputTag('TriggerResults', '', process.name_())
    if type(weight_src) == cms.InputTag:
        weight_src = weight_src.moduleLabel
    if weight_src:
        weight_obj = getattr(process, weight_src)
        process.epSimpleTriggerEfficiency = cms.EndPath(weight_obj * process.SimpleTriggerEfficiency)
    else:
        process.epSimpleTriggerEfficiency = cms.EndPath(process.SimpleTriggerEfficiency)

process.SimpleTriggerEfficiency.prescale_paths  = cms.vstring()  #*prescales.prescales.keys()),
process.SimpleTriggerEfficiency.prescale_values = cms.vuint32()  #*[o for l,h,o in prescales.prescales.itervalues()]),

process.SimpleTriggerEfficiencyMu      = process.SimpleTriggerEfficiency.clone()
process.SimpleTriggerEfficiencyMuInAcc = process.SimpleTriggerEfficiency.clone()
process.SimpleTriggerEfficiencyEl      = process.SimpleTriggerEfficiency.clone()
process.SimpleTriggerEfficiencyElInAcc = process.SimpleTriggerEfficiency.clone()

process.p1 = cms.Path(process.SimpleTriggerEfficiency)
#process.p2 = cms.Path(process.genMus      * process.genMuCount      * process.SimpleTriggerEfficiencyMu)
#process.p3 = cms.Path(process.genMusInAcc * process.genMuInAccCount * process.SimpleTriggerEfficiencyMuInAcc)
#process.p4 = cms.Path(process.genEls      * process.genElCount      * process.SimpleTriggerEfficiencyEl)
#process.p5 = cms.Path(process.genElsInAcc * process.genElInAccCount * process.SimpleTriggerEfficiencyElInAcc)
