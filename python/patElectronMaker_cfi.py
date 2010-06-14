import FWCore.ParameterSet.Config as cms

patElectronMaker = cms.EDProducer("PATElectronMaker",
	aliasPrefix = cms.untracked.string("els_pat"),
    # pat jet collection
    patElectronsInputTag  = cms.InputTag("selectedPatElectrons"),
    recoElectronsInputTag = cms.InputTag("gsfElectrons")
)


