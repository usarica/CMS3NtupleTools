import FWCore.ParameterSet.Config as cms

hltMaker = cms.EDProducer("HLTMaker",
    #processName = cms.untracked.string("HLT"),
    # HLT is the default process name, but it cannot be
    # guaranteed to be there. Use an empty string to
    # automatically choose the process name responsible
    # for the latest HLT products
    processName = cms.untracked.string(""),
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
