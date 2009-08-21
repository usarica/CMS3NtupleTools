import FWCore.ParameterSet.Config as cms

eventMaker = cms.EDFilter("EventMaker",
    exclusiveCrossSection = cms.untracked.double(0.0),
    kfactor = cms.untracked.double(1.0),
    inclusiveCrossSection = cms.untracked.double(0.0),
    haveL1TriggerInfo = cms.untracked.bool(True),
    haveHLTriggerInfo = cms.untracked.bool(True),
    datasetName = cms.string("undefined"),
    CMS2tag = cms.string("undefined")
)


