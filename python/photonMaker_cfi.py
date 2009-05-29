import FWCore.ParameterSet.Config as cms

photonMaker = cms.EDFilter("PhotonMaker",
    # Photon collection
    photonsInputTag = cms.InputTag("photons"),
    # Isolation
    ecalIsoTag = cms.InputTag("gamIsoFromDepsEcalFromHitsCMS2"),
    hcalIsoTag = cms.InputTag("gamIsoFromDepsHcalFromTowersCMS2"),
    tkIsoTag = cms.InputTag("gamIsoFromDepsTkCMS2")
)

