import FWCore.ParameterSet.Config as cms

jetMaker = cms.EDFilter("JetMaker",
           uncorJetsInputTag     = cms.InputTag("prunedUncorrectedCMS2Jets"),
           runningOnReco         = cms.untracked.bool(True),
           jetIDInputTag = cms.PSet(
                 hbheRecHitsColl = cms.InputTag("hbhereco"),
                 hoRecHitsColl   = cms.InputTag("horeco"),
                 hfRecHitsColl   = cms.InputTag("hfreco"),
                 ebRecHitsColl   = cms.InputTag("ecalRecHit", "EcalRecHitsEB"),
                 eeRecHitsColl   = cms.InputTag("ecalRecHit", "EcalRecHitsEE")
    )
)


