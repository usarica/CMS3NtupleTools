import FWCore.ParameterSet.Config as cms

electronMaker = cms.EDFilter("ElectronMaker",
    # Electron collection
    electronsInputTag = cms.InputTag("uniqueElectrons"),
    # Beamspot
    beamSpotInputTag = cms.InputTag("beamSpotMaker"),
    # Isolation
    ecalIsoTag = cms.InputTag("eleIsoFromDepsEcalFromHits"),
    hcalIsoTag = cms.InputTag("eleIsoFromDepsHcalFromHits"),
    tkIsoTag = cms.InputTag("eleIsoFromDepsTk")

)

