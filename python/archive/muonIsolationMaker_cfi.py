import FWCore.ParameterSet.Config as cms

muonIsolationMaker = cms.EDProducer("MuonIsolationMaker",
	
  aliasPrefix      = cms.untracked.string("mus"),
  muonsInputTag    = cms.InputTag("muons"        ),
  pfNoPileUpInputTag_ = cms.InputTag("pfNoPileUp"),
  cms2muonsInputTag = cms.InputTag("muonMaker","musp4")                           
)


