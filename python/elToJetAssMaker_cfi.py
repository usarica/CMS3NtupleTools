import FWCore.ParameterSet.Config as cms

elToJetAssMaker = cms.EDProducer("ElToJetAssMaker",
	aliasPrefix = cms.untracked.string("els"),
    # min DR
    minDR = cms.double(0.4),
    elsInputTag = cms.InputTag("electronMaker", "elsp4"),
    jetsInputTag = cms.InputTag("jetMaker", "jetsp4")                               
)
