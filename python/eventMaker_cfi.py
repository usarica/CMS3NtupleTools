import FWCore.ParameterSet.Config as cms

eventMaker = cms.EDFilter("EventMaker",
	aliasPrefix = cms.untracked.string("evt"),
    datasetName       = cms.string("undefined"),
    CMS2tag           = cms.string("V02-00-05"),
    dcsTag            = cms.InputTag("scalersRawToDigi")
)


