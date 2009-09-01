import FWCore.ParameterSet.Config as cms

photonMaker = cms.EDFilter("PhotonMaker",
    # Photon collection
    photonsInputTag = cms.InputTag("photons"),
)

