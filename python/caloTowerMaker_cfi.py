import FWCore.ParameterSet.Config as cms

caloTowerMaker = cms.EDFilter("CaloTowerMaker",
	aliasPrefix = cms.untracked.string("twrs"),
   primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
   caloTowersInputTag    = cms.InputTag("towerMaker"), #Default -- CLEANED
   #caloTowersInputTag    = cms.InputTag("cms2calotowermaker"), #from caloTowerSequence_cfi -- uncleaned
   ecalRecHitsInputTag_EE = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
   ecalRecHitsInputTag_EB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
   ecalDigiProducerEB = cms.InputTag("ecalDigis:ebDigis"),
   ecalDigiProducerEE = cms.InputTag("ecalDigis:eeDigis"),
   threshEt       = cms.double(5.), #gev, for reading out flags, time, etc
   spikeR4Thresh  = cms.double(0.05), #spike if s4/s1 < this and et > next
   spikeEtThresh  = cms.double(5.), #gev
   spikeEtaMax    = cms.double(1.4442), #exclude edge of barrel
)

