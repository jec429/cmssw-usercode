import FWCore.ParameterSet.Config as cms

ana = cms.EDAnalyzer('WR_Analyzer',
                    gen_src = cms.InputTag('prunedGenParticles'),
                    vertex_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
                    gen_jet_src = cms.InputTag('slimmedGenJets'),
                    electron_src = cms.InputTag('slimmedElectrons'),
                    electronVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-veto"),
                    electronLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-loose"),
                    electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-medium"),
                    electronTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-tight"),
                    electronHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51-miniAOD"),                            
                    muon_src = cms.InputTag('slimmedMuons'),
                    jet_src = cms.InputTag('slimmedJets'),
                    jet_pt_min = cms.double(40),
                    jet_eta_max = cms.double(2.5),
                    dilepton_mass_cut = cms.double(0.0),
                    lljj_mass_cut = cms.double(0.0)
                    )
