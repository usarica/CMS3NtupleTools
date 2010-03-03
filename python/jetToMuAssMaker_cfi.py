import FWCore.ParameterSet.Config as cms

jetToMuAssMaker = cms.EDFilter("JetToMuAssMaker",
    # min DR
    minDR = cms.double(0.4),
    jetsInputTag = cms.InputTag("jetMaker", "jetsp4"),
    musInputTag  = cms.InputTag("muonMaker","musp4")
)


