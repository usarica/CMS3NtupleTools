import FWCore.ParameterSet.Config as cms

photonMaker = cms.EDFilter("PhotonMaker",
	aliasPrefix = cms.untracked.string("photons"),
    minEt       = cms.double(10.), #gev, min to keep
    # Photon collection
    photonsInputTag = cms.InputTag("photons"),
    ecalRecHitsInputTag_EE = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    ecalRecHitsInputTag_EB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    cms2scsseeddetidInputTag = cms.InputTag("scMaker"),
)

