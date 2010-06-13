import FWCore.ParameterSet.Config as cms

pfElToElAssMaker = cms.EDFilter("PFElToElAssMaker",
                                aliasPrefix = cms.untracked.string("pfels"),                                
                                minDR = cms.double(0.1),
                                elsInputTag = cms.InputTag("electronMaker"),
                                pfelsInputTag = cms.InputTag("pfElectronMaker")
)


