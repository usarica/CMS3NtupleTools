import FWCore.ParameterSet.Config as cms

eventMaker = cms.EDFilter("EventMaker",
    exclusiveCrossSection = cms.untracked.double(0.0),
    kfactor = cms.untracked.double(1.0),
    inclusiveCrossSection = cms.untracked.double(0.0),
    haveTriggerInfo = cms.untracked.bool(True)
)


