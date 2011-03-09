import FWCore.ParameterSet.Config as cms

fastJetMaker = cms.EDProducer("FastJetMaker",
                              aliasPrefix = cms.untracked.string("evt"),
                              rho_tag     = cms.InputTag("kt6PFJets","rho")
                              )
