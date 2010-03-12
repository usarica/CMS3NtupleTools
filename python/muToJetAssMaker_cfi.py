import FWCore.ParameterSet.Config as cms

muToJetAssMaker = cms.EDFilter("MuToJetAssMaker",
	aliasPrefix = cms.untracked.string("mus"),
    # min DR
    minDR = cms.double(0.4),
    musInputTag  = cms.InputTag("muonMaker"),
    jetsInputTag = cms.InputTag("jetMaker")                                
)


