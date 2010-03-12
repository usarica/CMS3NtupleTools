import FWCore.ParameterSet.Config as cms

elToJetAssMaker = cms.EDFilter("ElToJetAssMaker",
	aliasPrefix = cms.untracked.string("els"),
    # min DR
    minDR = cms.double(0.4),
    elsInputTag = cms.InputTag("electronMaker"),
    jetsInputTag = cms.InputTag("jetMaker")                               
)
