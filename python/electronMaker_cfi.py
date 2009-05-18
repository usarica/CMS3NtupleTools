import FWCore.ParameterSet.Config as cms

electronMaker = cms.EDFilter("ElectronMaker",
    # Jet collection
    electronsInputTag = cms.InputTag("pixelMatchGsfElectrons"),
    # Beamspot
    beamSpotInputTag = cms.InputTag("beamSpotMaker"),
    # Isolation
    ecalIsoTag = cms.InputTag("eleIsoFromDepsEcalFromHits"),
    hcalIsoTag = cms.InputTag("eleIsoFromDepsHcalFromHits"),
    tkIsoTag = cms.InputTag("eleIsoFromDepsTk")

)

