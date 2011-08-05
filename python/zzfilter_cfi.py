import FWCore.ParameterSet.Config as cms

fourLeptons = cms.EDFilter("ZZFilter",
    muons           = cms.InputTag("muons"),
    electrons       = cms.InputTag("gsfElectrons"),
    minMuPt         = cms.double(5.),
    minElePt        = cms.double(5.),
    minLeadingLepPt = cms.double(20.),                           
    minNLeptons     = cms.int32(4),
    prescale        = cms.int32(1)
)                                              


