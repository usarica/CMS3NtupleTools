import FWCore.ParameterSet.Config as cms

caloTowerHFMaker = cms.EDFilter("CaloTowerHFMaker",
   aliasPrefix = cms.untracked.string("twrs"),
   caloTowersInputTag    = cms.InputTag("towerMaker"),
   cms2TowersInputTag    = cms.InputTag("caloTowerMaker"),
   hfReflaggedHitsInputTag  = cms.InputTag("hfrecoReflagged"),
   threshHF     = cms.double(5.),							  
)

