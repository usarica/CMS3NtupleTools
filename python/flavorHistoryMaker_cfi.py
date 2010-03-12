import FWCore.ParameterSet.Config as cms

flavorHistoryMaker = cms.EDFilter("FlavorHistoryMaker",                                  
	aliasPrefix = cms.untracked.string("genps"),
    flavorHistoryFilterTag  = cms.InputTag("flavorHistoryFilter")
)

