import FWCore.ParameterSet.Config as cms

ele10mu10IsoIdMET = cms.EDFilter("WWFilter",
    muons           = cms.InputTag("muons"),
    electrons       = cms.InputTag("gsfElectrons"),
    mets            = cms.VInputTag( cms.InputTag("pfMet"), cms.InputTag("tcMet") ),
    minMuPt         = cms.double(10.),
    minElePt        = cms.double(10.),
    minMass         = cms.double(5.),
    minMET          = cms.double(20.),
    applyMuId       = cms.bool(True),
    applyMuIso      = cms.bool(True),
    applyEleId      = cms.bool(True),
    applyEleIso     = cms.bool(True),
    prescale        = cms.int32(1)
)                                              
                                              
ele10mu10IsoId = cms.EDFilter("WWFilter",
    muons           = cms.InputTag("muons"),
    electrons       = cms.InputTag("gsfElectrons"),
    mets            = cms.VInputTag( cms.InputTag("pfMet"), cms.InputTag("tcMet") ),
    minMuPt         = cms.double(10.),
    minElePt        = cms.double(10.),
    minMass         = cms.double(5.),
    minMET          = cms.double(-999.),
    applyMuId       = cms.bool(True),
    applyMuIso      = cms.bool(True),
    applyEleId      = cms.bool(True),
    applyEleIso     = cms.bool(True),
    prescale        = cms.int32(100)
)                                              


