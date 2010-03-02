import FWCore.ParameterSet.Config as cms

patElectronMaker = cms.EDFilter("PATElectronMaker",
    # pat jet collection
    patElectronsInputTag  = cms.InputTag("selectedPatElectrons"),
    recoElectronsInputTag = cms.InputTag("gsfElectrons")
)


