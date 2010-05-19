import FWCore.ParameterSet.Config as cms


ecalRecHitMaker = cms.EDFilter("EcalRecHitMaker",
                               ecalEBRecHitInputTag = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
                               ecalEERecHitInputTag = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
                               minEt                = cms.double(5.),
                               AliasPrefix          = cms.string("ecalrhit"),
                               )


