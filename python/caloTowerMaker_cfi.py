import FWCore.ParameterSet.Config as cms

caloTowerMaker = cms.EDFilter("CaloTowerMaker",
   primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
   caloTowersInputTag    = cms.InputTag("towerMaker"),
   ecalRecHitsInputTag_EE = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
   ecalRecHitsInputTag_EB = cms.InputTag("ecalRecHit","EcalRecHitsEB")

#    scEtMin = cms.double(10.0)

)

