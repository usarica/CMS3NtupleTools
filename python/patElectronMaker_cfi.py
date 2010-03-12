import FWCore.ParameterSet.Config as cms

patElectronMaker = cms.EDFilter("PATElectronMaker",
	aliasPrefix = cms.untracked.string("els_pat"),
    # pat jet collection
    patElectronsInputTag  = cms.InputTag("selectedPatElectrons"),
    recoElectronsInputTag = cms.InputTag("gsfElectrons")
)


