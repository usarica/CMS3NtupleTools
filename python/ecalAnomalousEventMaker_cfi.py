import FWCore.ParameterSet.Config as cms

ecalAnomalousEventMaker = cms.EDProducer(
  "EcalAnomalousEventMaker",
  aliasPrefix = cms.untracked.string("evt")
)


