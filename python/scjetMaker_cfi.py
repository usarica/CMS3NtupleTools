import FWCore.ParameterSet.Config as cms

scjetMaker = cms.EDFilter("JetMaker",
           uncorJetsInputTag     = cms.InputTag("prunedUncorrectedCMS2scJets"),
           runningOnReco         = cms.untracked.bool(True),
           jetIDInputTag = cms.PSet(
                 hbheRecHitsColl = cms.InputTag("hbhereco"),
                 hoRecHitsColl   = cms.InputTag("horeco"),
                 hfRecHitsColl   = cms.InputTag("hfreco"),
                 ebRecHitsColl   = cms.InputTag("ecalRecHit", "EcalRecHitsEB"),
                 eeRecHitsColl   = cms.InputTag("ecalRecHit", "EcalRecHitsEE")
                 ),
           L2L3JetCorrectorName = cms.string("L2L3JetCorrectorSC5Calo"),
           AliasPrefix          = cms.string("scjets")
)


