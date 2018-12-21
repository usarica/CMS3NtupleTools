import FWCore.ParameterSet.Config as cms
# import CMS3.NtupleMaker.configProcessName as configProcessName

metFilterMaker = cms.EDProducer(
  "MetFilterMaker",
  aliasPrefix      = cms.untracked.string("filt"),
  filtersInputTag  = cms.InputTag("TriggerResults"),
  # processName      = cms.untracked.string(configProcessName.name)
  doEcalFilterUpdate = cms.bool(True),
)
