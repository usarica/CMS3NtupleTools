import FWCore.ParameterSet.Config as cms

hltMaker = cms.EDProducer("HLTMaker",
   aliasprefix = cms.untracked.string("hlt"),
   processName = cms.untracked.string("HLT"),
   prunedTriggerNames = cms.untracked.vstring(
      "HLT*Mu*",
      "HLT*Ele*",
      "HLT*EG*",
      "HLT*Photon*",
      "HLT*Jet*",
      "HLT*MET*",
   ),
   triggerPrescaleInputTag = cms.untracked.string("patTrigger"),
   #fillTriggerObjects = cms.untracked.bool(False),
   #triggerObjectsName = cms.untracked.string("selectedPatTrigger"),
)

