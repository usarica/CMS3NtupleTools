import FWCore.ParameterSet.Config as cms

eventMaker = cms.EDProducer(
  "EventMaker",
  aliasPrefix = cms.untracked.string("evt"),
  datasetName = cms.string("undefined"),
  CMS3tag     = cms.string("V02-00-05"),
  dcsTag      = cms.InputTag("scalersRawToDigi"),
  scalersTag  = cms.InputTag("scalersRawToDigi"),
  isData      = cms.bool(False)
)


