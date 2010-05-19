import FWCore.ParameterSet.Config as cms

#There are three different HLT trigger menus in the Spring10 samples,
#REDIGI, HLT and HLT8e29. The appropriate one is REDIGI, and this will be
#the default although we will run all 3 paths


hlt8e29Maker = cms.EDProducer("HLTMaker",
    processName = cms.untracked.string("HLT8E29"),
    aliasPrefix = cms.untracked.string("hlt8e29"),                       
    fillTriggerObjects = cms.untracked.bool(True),
    prunedTriggerNames = cms.untracked.vstring(
        "HLT*Mu*",
        "HLT*Ele*",
        "HLT*EG*",
        "HLT*Photon*",
        "HLT*Jet*",
        "HLT*MET*",
    )
)

hlt1e31Maker = cms.EDProducer("HLTMaker",
    processName = cms.untracked.string("HLT"),
    aliasPrefix = cms.untracked.string("hlt1e31"),
    fillTriggerObjects = cms.untracked.bool(True),
    prunedTriggerNames = cms.untracked.vstring(
        "HLT*Mu*",
        "HLT*Ele*",
        "HLT*EG*",
        "HLT*Photon*",
        "HLT*Jet*",
        "HLT*MET*",
    )
)

hltMaker = cms.EDProducer("HLTMaker",
    processName = cms.untracked.string("REDIGI"),
    aliasPrefix = cms.untracked.string("hlt"),                       
    fillTriggerObjects = cms.untracked.bool(True),
    prunedTriggerNames = cms.untracked.vstring(
        "HLT*Mu*",
        "HLT*Ele*",
        "HLT*EG*",
        "HLT*Photon*",
        "HLT*Jet*",
        "HLT*MET*",
    )
)
