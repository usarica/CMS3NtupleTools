import FWCore.ParameterSet.Config as cms

ecalAnomalousEventMaker = cms.EDProducer(
  "EcalAnomalousEventMaker",
  aliasPrefix          = cms.untracked.string("evt"),
  JetCorrectionService = cms.string('ak5CaloL1L2L3Residual') ## or 'ak5PFL1L2L3Residual' or 'ak5JPTL1L2L3Residual'
)


