import FWCore.ParameterSet.Config as cms

metFilterMaker = cms.EDProducer(
  "MetFilterMaker",
  aliasPrefix      = cms.untracked.string("filt"),
  filtersInputTag  = cms.InputTag("TriggerResults"),
  processName      = cms.untracked.string("HLT" )

)
