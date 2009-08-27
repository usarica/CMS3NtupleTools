import FWCore.ParameterSet.Config as cms

elCaloIsoMaker = cms.EDProducer("ElCaloIsoMaker",
  electronsInputTag = cms.InputTag("gsfElectrons"),
  basicClusterInputTag = cms.InputTag("egammaBasicClusterMerger"),
  caloTowersInputTag = cms.InputTag("towerMaker"),
  minDR = cms.double(0.05),
  maxDR = cms.double(0.4),
  minDEta = cms.double(0.01)
)
