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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring('file:myFile.root'),
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
                        
process.source.fileNames = file_list(dset+'step4_*.root',True)
#process.source.fileNames = file_list('/WRToNuEToEEJJ_MW-2000_MNu-1000_TuneCUETP8M1_13TeV-pythia8/RunIIFall14GS-MCRUN2_71_V1-v1/GEN-SIM',False)
outfile = 'comparison.root'

if ttbar:
    process.source.fileNames = file_list(Samples.ttbar_samples[0].dataset,False)
    outfile = 'comparison_ttbar.root'
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

if dyjets:
    process.source.fileNames = file_list(Samples.dyjets_samples[0].dataset,False)
    outfile = 'comparison_dyjets.root'
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

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
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_cff']
#Add them to the VID producer
#for idmod in my_id_modules:
    #setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
process.load('RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff')
process.heeptest2 = process.heepElectronID_HEEPV51_miniAOD.clone(idName = cms.string("heeptest2"))
process.heeptest2.cutFlow[0].isIgnored = True
process.heepElectronID_HEEPV51_miniAOD.cutFlow[0].minPt = 20.0

setupVIDSelection(process.egmGsfElectronIDs,process.heepElectronID_HEEPV51_miniAOD)
setupVIDSelection(process.egmGsfElectronIDs,process.heeptest2)

###################################################################################################
###################################################################################################

process.ana = cms.EDAnalyzer('Comparison',
                             electronHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51-miniAOD")
                             #electronHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heeptest2")
                             )

process.p = cms.Path(process.egmGsfElectronIDSequence * process.ana)
