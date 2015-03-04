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
                            fileNames=cms.untracked.vstring('file:EXO-Phys14DR-00009.root'),
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
#process.source.fileNames = file_list('/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',False)
#process.source.fileNames = file_list('/eos/uscms/store/user/jchaves/Nstep_EE_2000_reco/EXO-Phys14DR-00009_*.root',True)
process.source.fileNames = file_list('/eos/uscms/store/user/jchaves/Nstep_MUMU_2000_reco/EXO-Phys14DR-00009_*.root',True)
outfile = 'out_EXO.root'

if ttbar:
    process.source.fileNames = file_list(Samples.ttbar_samples[0].dataset,False)
    outfile = 'out_ttbar.root'
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

if dyjets:
    process.source.fileNames = file_list(Samples.dyjets_samples[0].dataset,False)
    outfile = 'out_dyjets.root'
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.TFileService = cms.Service('TFileService', fileName = cms.string(outfile))

    
###################################################################################################
###################################################################################################
##########################       Electron ID                   ####################################
###################################################################################################
###################################################################################################

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
#    for PHYS14 scenario PU4bx50 : global tag is ???
#    for PHYS14 scenario PU20bx25: global tag is PHYS14_25_V1
#  as a rule, find the global tag in the DAS under the Configs for given dataset
process.GlobalTag.globaltag = 'PHYS14_25_V1'

# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
# overwrite a default parameter: for miniAOD, the collection name is a slimmed one
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')

from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)

# Define which IDs we want to produce
# Each of these two example IDs contains all four standard 
# cut-based ID working points (only two WP of the PU20bx25 are actually used here).
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_cff']
#Add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
process.load('RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff')
process.heeptest2 = process.heepElectronID_HEEPV51_miniAOD.clone(idName = cms.string("heeptest2"))
process.heeptest2.cutFlow[0].isIgnored = True
process.heeptest2.cutFlow[1].isIgnored = True
process.heeptest2.cutFlow[2].isIgnored = True
process.heeptest2.cutFlow[3].isIgnored = True
process.heeptest2.cutFlow[4].isIgnored = True
process.heeptest2.cutFlow[5].isIgnored = True
process.heeptest2.cutFlow[6].isIgnored = True
process.heeptest2.cutFlow[7].isIgnored = True
process.heeptest2.cutFlow[8].isIgnored = True
process.heeptest2.cutFlow[9].isIgnored = True
process.heeptest2.cutFlow[10].isIgnored = True
process.heeptest2.cutFlow[0].minPt = 100.0

setupVIDSelection(process.egmGsfElectronIDs,process.heepElectronID_HEEPV51_miniAOD)
setupVIDSelection(process.egmGsfElectronIDs,process.heeptest2)

###################################################################################################
###################################################################################################

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.triggerFilter = hltHighLevel.clone()
process.triggerFilter.HLTPaths = ['HLT_Photon155_v*','HLT_Mu40_v*']
process.triggerFilter.andOr = True # = OR
#for name, path in process.paths.items():
 #   if not name.startswith('eventCleaning'):
  #      path.insert(0, process.triggerFilter)
process.ptrig = cms.Path(process.triggerFilter)

process.load('JChaves.WR_Analyzer.eventSelector_cfi')
process.load('JChaves.WR_Analyzer.analyzer_cfi')

process.fltr.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('ptrig'))

#process.ana2 = process.ana.clone(electronHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heeptest2"))
process.ana2 = process.ana.clone(dilepton_mass_cut = cms.double(200.0))
process.ana3 = process.ana.clone(dilepton_mass_cut = cms.double(200.0),
                                   lljj_mass_cut = cms.double(600.0))
process.ana4 = process.ana.clone(dilepton_mass_cut = cms.double(200.0),
                                   lljj_mass_cut = cms.double(600.0),
                                   first_l_pt_cut = cms.double(100.0))

process.p = cms.Path(process.egmGsfElectronIDSequence * process.fltr * process.ana)
process.p2 = cms.Path(process.egmGsfElectronIDSequence * process.fltr * process.ana2)
process.p3 = cms.Path(process.egmGsfElectronIDSequence * process.fltr * process.ana3)
process.p4 = cms.Path(process.egmGsfElectronIDSequence * process.fltr * process.ana4)

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
        for x in ['HT-100to200','HT-200to400','HT-400to600','HT-600toInf']:
        #for x in ['M-200To400','M-400To800','M-800To1400','M-1400To2300','M-3500To4500','M-4500To6000','M-6000To7500','M-7500To8500',]:
            #crab2_submit('dyjets',-1,x)
            crab3_submit('dyjets',-1,x)
