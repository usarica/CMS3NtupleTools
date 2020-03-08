import FWCore.ParameterSet.Config as cms

hltMaker = cms.EDProducer("HLTMaker",
   aliasprefix = cms.untracked.string("hlt"),
   processName = cms.untracked.string("HLT"),
   prunedTriggerNames = cms.untracked.vstring(),
   triggerPrescaleInputTag = cms.untracked.string("patTrigger"),
   recordFilteredTrigObjects = cms.bool(False),
   triggerObjectsName = cms.untracked.string("slimmedPatTrigger"),
   )


