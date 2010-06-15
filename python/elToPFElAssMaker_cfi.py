import FWCore.ParameterSet.Config as cms

elToPFElAssMaker = cms.EDProducer("ElToPFElAssMaker",
                                aliasPrefix = cms.untracked.string("els"),                                
                                minDR = cms.double(0.1),
                                elsInputTag = cms.InputTag("electronMaker"),
                                pfelsInputTag = cms.InputTag("pfElectronMaker")
)


