import FWCore.ParameterSet.Config as cms

photonMaker = cms.EDFilter("PhotonMaker",
	aliasPrefix = cms.untracked.string("photons"),
    # Photon collection
    photonsInputTag = cms.InputTag("photons"),
)

