import FWCore.ParameterSet.Config as cms

fastJetMaker = cms.EDProducer("FastJetMaker",
                              aliasPrefix = cms.untracked.string("evt"),
                              rhoJEC_tag     = cms.InputTag("kt6PFJetsDeterministicJEC","rho"),
                              rhoIso_tag     = cms.InputTag("kt6PFJetsDeterministicIso","rho")
                              )

wwRhoDefault = cms.EDProducer("EnergyDensityMaker",
                              input = cms.InputTag("kt6PFJetsForRhoComputationDefault","rho"),
                              alias = cms.untracked.string("evt_ww_rho"))
                              
wwRhoVoronoi = cms.EDProducer("EnergyDensityMaker",
                              input = cms.InputTag("kt6PFJetsForRhoComputationVoronoi","rho"),
                              alias = cms.untracked.string("evt_ww_rho_vor"))

wwRhoRandom = cms.EDProducer("EnergyDensityMaker",
                             input = cms.InputTag("kt6PFJetsForRhoComputationRandom","rho"),
                             alias = cms.untracked.string("evt_ww_rho_rnd"))
