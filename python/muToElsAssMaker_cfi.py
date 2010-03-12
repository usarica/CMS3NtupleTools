import FWCore.ParameterSet.Config as cms

muToElsAssMaker = cms.EDFilter("MuToElsAssMaker",
	aliasPrefix = cms.untracked.string("mus"),
    # min DR
    minDR = cms.double(0.1),
    musInputTag = cms.InputTag("muonMaker"),
    elsInputTag = cms.InputTag("electronMaker")
)


