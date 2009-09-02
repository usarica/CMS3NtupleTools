import FWCore.ParameterSet.Config as cms

hltMaker8e29 = cms.EDProducer("HLTMaker",
    processName = cms.untracked.string("HLT8E29"),
    fillTriggerObjects = cms.untracked.bool(True),
    prunedTriggerNames = cms.untracked.vstring("HLT_Mu9")
)

hltMaker1e31 = cms.EDProducer("HLTMaker",
    processName = cms.untracked.string("HLT"),
    fillTriggerObjects = cms.untracked.bool(True),
    prunedTriggerNames = cms.untracked.vstring("HLT_Mu11")
)
