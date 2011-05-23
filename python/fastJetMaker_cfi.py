import FWCore.ParameterSet.Config as cms

fastJetMaker = cms.EDProducer("FastJetMaker",
                              aliasPrefix = cms.untracked.string("evt"),
                              rhoJEC_tag     = cms.InputTag("kt6PFJetsDeterministicJEC","rho"),
                              rhoIso_tag     = cms.InputTag("kt6PFJetsDeterministicIso","rho")
                              )
