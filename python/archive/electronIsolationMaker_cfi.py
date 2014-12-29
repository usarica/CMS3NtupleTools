import FWCore.ParameterSet.Config as cms

electronIsolationMaker = cms.EDProducer("ElectronIsolationMaker",	
                                        aliasPrefix      = cms.untracked.string("els"),
                                        gsfElectronInputTag    = cms.InputTag("gedGsfElectrons"),
                                        pfNoPileUpInputTag_ = cms.InputTag("pfNoPileUp"),
                                        cms2electronInputTag = cms.InputTag("electronMaker","elsp4"),
                                        vertexInputTag_ = cms.InputTag("offlinePrimaryVertices")
)


