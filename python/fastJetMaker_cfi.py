import FWCore.ParameterSet.Config as cms

fastJetMaker = cms.EDProducer("FastJetMaker",
                              aliasPrefix = cms.untracked.string("evt"),
                              rhoJEC_tag     = cms.InputTag("kt6PFJetsDeterministicJEC","rho"),
                              rhoIso_tag     = cms.InputTag("kt6PFJetsDeterministicIso","rho")
                              )

wwRhoDefaultMaker = cms.EDProducer("EnergyDensityMaker",
                              input = cms.InputTag("kt6PFJets","rho"),
                              alias = cms.untracked.string("evt_ww_rho"))
                              
wwRhoVoronoiMaker = cms.EDProducer("EnergyDensityMaker",
                              input = cms.InputTag("kt6PFJetsForRhoComputationVoronoi","rho"),
                              alias = cms.untracked.string("evt_ww_rho_vor"))
                              
wwRhoActiveAreaMaker = cms.EDProducer("EnergyDensityMaker",
                              input = cms.InputTag("kt6PFJetsForRhoComputationActiveArea","rho"),
                              alias = cms.untracked.string("evt_ww_rho_act"))

wwRhoRandomMaker = cms.EDProducer("EnergyDensityMaker",
                             input = cms.InputTag("kt6PFJetsForRhoComputationRandom","rho"),
                             alias = cms.untracked.string("evt_ww_rho_rnd"))

fixedGridRhoAllMaker = cms.EDProducer("EnergyDensityMaker",
                                      input = cms.InputTag("fixedGridRhoAll"),
                                      alias = cms.untracked.string("evt_fixgrid_all_rho"))

fixedGridRhoFastJetAllMaker = cms.EDProducer("EnergyDensityMaker",
                                             input = cms.InputTag("fixedGridRhoFastjetAll"),
                                             alias = cms.untracked.string("evt_fixgridfastjet_all_rho"))

kt6CaloJetsRhoMaker = cms.EDProducer("EnergyDensityMaker",
                                     input = cms.InputTag("kt6CaloJets","rho"),
                                     alias = cms.untracked.string("evt_kt6calo_rho"))

kt6CaloJetsCentralRhoMaker = cms.EDProducer("EnergyDensityMaker",
                                            input = cms.InputTag("kt6CaloJetsCentral","rho"),
                                            alias = cms.untracked.string("evt_kt6calo_central_rho"))

kt6PFJetsCentralChargedPileUpRhoMaker = cms.EDProducer("EnergyDensityMaker",
                                                       input = cms.InputTag("kt6PFJetsCentralChargedPileUp","rho"),
                                                       alias = cms.untracked.string("evt_kt6pf_ctrChargedPU_rho"))

kt6PFJetsCentralNeutralRhoMaker = cms.EDProducer("EnergyDensityMaker",
                                                 input = cms.InputTag("kt6PFJetsCentralNeutral","rho"),
                                                 alias = cms.untracked.string("evt_kt6pf_ctrNeutral_rho"))

kt6PFJetsCentralNeutralTightRhoMaker = cms.EDProducer("EnergyDensityMaker",
                                                      input = cms.InputTag("kt6PFJetsCentralNeutralTight","rho"),
                                                      alias = cms.untracked.string("evt_kt6pf_ctrNeutralTight_rho"))


kt6PFJetsForEGIsolationRhoMaker = cms.EDProducer("EnergyDensityMaker",
                                                 input = cms.InputTag("kt6PFJetsForEGIsolation","rho"),
                                                 alias = cms.untracked.string("evt_kt6pf_foregiso_rho"))
