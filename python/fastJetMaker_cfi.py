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

wwRhoRandomMaker = cms.EDProducer("EnergyDensityMaker",
                             input = cms.InputTag("kt6PFJetsForRhoComputationRandom","rho"),
                             alias = cms.untracked.string("evt_ww_rho_rnd"))
