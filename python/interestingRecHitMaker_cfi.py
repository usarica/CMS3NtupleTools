import FWCore.ParameterSet.Config as cms

interestingRecHitMaker = cms.EDFilter("InterestingRecHitMaker",
aliasPrefix = cms.untracked.string("interestingrechit"),
scDetIdCMS2 = cms.InputTag("scMaker", "scsdetIdSeed"),
ecalRecHitsInputTag_EE = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
ecalRecHitsInputTag_EB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
ecalDigiProducerEB = cms.InputTag("ecalDigis:ebDigis"),
ecalDigiProducerEE = cms.InputTag("ecalDigis:eeDigis"),
threshEt   = cms.double(5.), #gev, for reading out flags, time, etc
spikeR4Thresh  = cms.double(0.05), #spike if s4/s1 < this and et > next
spikeEtThresh  = cms.double(5.), #gev
spikeEtaMax= cms.double(1.4442), #exclude edge of barrel
)

