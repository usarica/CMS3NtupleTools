import FWCore.ParameterSet.Config as cms

electronIsolationMaker = cms.EDProducer("ElectronIsolationMaker",
	
  aliasPrefix      = cms.untracked.string("els"),
  gsfElectronInputTag    = cms.InputTag("gsfElectrons"),
  pfNoPileUpInputTag_ = cms.InputTag("pfNoPileUp"),
  cms2electronInputTag = cms.InputTag("electronMaker","elsp4")                           
)


