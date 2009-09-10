import FWCore.ParameterSet.Config as cms

eventMaker = cms.EDFilter("EventMaker",
    datasetName       = cms.string("undefined"),
    CMS2tag           = cms.string("V02-00-01")
)


