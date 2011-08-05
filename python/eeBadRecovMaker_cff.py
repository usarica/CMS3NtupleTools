import FWCore.ParameterSet.Config as cms

eeBadRecovMaker = cms.EDProducer("EEBadRecovMaker",

  ecalEBRecHitInputTag = cms.InputTag("ecalRecHit","reducedEcalRecHitsEB"),
  ecalEERecHitInputTag = cms.InputTag("ecalRecHit","reducedEcalRecHitsEE")
                               
)


