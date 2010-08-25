import FWCore.ParameterSet.Config as cms

#flavorHistoryMaker = cms.EDFilter("FlavorHistoryMaker",                                  
flavorHistoryMaker = cms.EDProducer("FlavorHistoryMaker",                                  
	aliasPrefix = cms.untracked.string("genps"),
    flavorHistoryFilterTag  = cms.InputTag("flavorHistoryFilter")
)

