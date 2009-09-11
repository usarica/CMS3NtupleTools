import FWCore.ParameterSet.Config as cms

caloTowerMaker = cms.EDFilter("CaloTowerMaker",
   primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
   caloTowersInputTag    = cms.InputTag("towerMaker")
#    scEtMin = cms.double(10.0)
)

