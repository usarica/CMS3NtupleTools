import FWCore.ParameterSet.Config as cms

patElectronMaker = cms.EDFilter("PATElectronMaker",
    # pat jet collection
    patElectronsInputTag  = cms.InputTag("selectedLayer1Electrons"),
    recoElectronsInputTag = cms.InputTag("gsfElectrons")
)


