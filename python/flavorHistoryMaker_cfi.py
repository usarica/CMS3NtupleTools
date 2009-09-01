import FWCore.ParameterSet.Config as cms

flavorHistoryMaker = cms.EDFilter("FlavorHistoryMaker",                                  
    flavorHistoryFilterTag  = cms.InputTag("flavorHistoryFilter")
)

