import FWCore.ParameterSet.Config as cms

muToElsAssMaker = cms.EDProducer("MuToElsAssMaker",
	aliasPrefix = cms.untracked.string("mus"),
    # min DR
    minDR = cms.double(0.1),
    musInputTag = cms.InputTag("muonMaker", "musp4"),
    elsInputTag = cms.InputTag("electronMaker", "elsp4")
)


