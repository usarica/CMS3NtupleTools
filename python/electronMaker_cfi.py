import FWCore.ParameterSet.Config as cms

electronMaker = cms.EDFilter("ElectronMaker",
    # Electron collection
    electronsInputTag = cms.InputTag("uniqueElectrons"),
    # Beamspot
    beamSpotInputTag = cms.InputTag("beamSpotMaker"),
    # Isolation
    ecalIsoTag = cms.InputTag("eleIsoFromDepsEcalFromHitsCMS2"),
    hcalIsoTag = cms.InputTag("eleIsoFromDepsHcalFromTowersCMS2"),
    tkIsoTag = cms.InputTag("eleIsoFromDepsTkCMS2")
)

