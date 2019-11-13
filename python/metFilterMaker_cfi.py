import FWCore.ParameterSet.Config as cms

metFilterMaker = cms.EDProducer(
  "METFilterMaker",
  aliasprefix      = cms.untracked.string("metfilt"),
  doEcalFilterUpdate = cms.bool(True),

  filtersInputTag  = cms.InputTag("TriggerResults"),
  # processName      = cms.untracked.string(configProcessName.name)
)
