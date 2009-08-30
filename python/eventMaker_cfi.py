import FWCore.ParameterSet.Config as cms

eventMaker = cms.EDFilter("EventMaker",
    haveL1TriggerInfo = cms.untracked.bool(True),
    haveHLTriggerInfo = cms.untracked.bool(True),
    datasetName       = cms.string("undefined"),
    CMS2tag           = cms.string("undefined")
)


