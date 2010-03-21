import FWCore.ParameterSet.Config as cms

spikeFilter = cms.EDFilter("TheFilter",
    # this syntax means: keep events that contain more than 0
    # entries in the collection caloTowerMaker:spikeEt -- spikeEt and spikeR4 should always have same n-entries
    filterExpressions = cms.VInputTag(cms.InputTag("caloTowerMaker","twrsspikeEt"))
)


