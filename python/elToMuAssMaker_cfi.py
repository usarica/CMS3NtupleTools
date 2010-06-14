import FWCore.ParameterSet.Config as cms

elToMuAssMaker = cms.EDProducer("ElToMuAssMaker",
	aliasPrefix = cms.untracked.string("els"),
    # min DR
    minDR = cms.double(0.1),
    elsInputTag = cms.InputTag("electronMaker"),
    musInputTag = cms.InputTag("muonMaker")
)


