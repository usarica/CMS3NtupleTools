#This is to rerun the "CaloTowersCreator" (see  RecoLocalCalo/ CaloTowersCreator/ python/ calotowermaker_cfi.py) WITHOUT cleaning of spikes

import FWCore.ParameterSet.Config as cms

from CMS2.NtupleMaker.caloTowerMaker_cfi import *
from RecoLocalCalo.CaloTowersCreator.calotowermaker_cfi import *
from Geometry.CaloEventSetup.CaloTowerConstituents_cfi import *
from RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi import *

cms2calotowermaker = calotowermaker.clone() 

#default is 1--excludes spikes which are kWeird==3. 4 should include them, but exclude kBad==4.
cms2calotowermaker.EcalAcceptSeverityLevel = cms.uint32(4)
#leave these off to avoid too many changes/complications
#cms2calotowermaker.UseEcalRecoveredHits = cms.bool(True)

#this is a temporary measure until the above does not exclude anything which would be below
caloTowerMakerUncleaned = cms.EDFilter("CaloTowerMaker",
   aliasPrefix = cms.untracked.string("twrsUncleaned"),
   primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
   #caloTowersInputTag    = cms.InputTag("towerMaker"), #default -- CLEANED
   caloTowersInputTag    = cms.InputTag("cms2calotowermaker"), #from caloTowerSequence_cfi -- UNCLEANED
   ecalRecHitsInputTag_EE = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
   ecalRecHitsInputTag_EB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
   ecalDigiProducerEB = cms.InputTag("ecalDigis:ebDigis"),
   ecalDigiProducerEE = cms.InputTag("ecalDigis:eeDigis"),
   threshEt       = cms.double(5.), #gev, for reading out flags, time, etc
   spikeR4Thresh  = cms.double(0.05), #spike if s4/s1 < this and et > next
   spikeEtThresh  = cms.double(5.), #gev
   spikeEtaMax    = cms.double(1.4442), #exclude edge of barrel
)


cms2CaloTowerSequence = cms.Sequence( cms2calotowermaker * caloTowerMaker * caloTowerMakerUncleaned )


