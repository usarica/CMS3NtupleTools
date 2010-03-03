import FWCore.ParameterSet.Config as cms

muToElsAssMaker = cms.EDFilter("MuToElsAssMaker",
    # min DR
    minDR = cms.double(0.1),
    musInputTag = cms.InputTag("muonMaker"),
    elsInputTag = cms.InputTag("electronMaker")
)


