import FWCore.ParameterSet.Config as cms

duplicateFilter = cms.EDFilter("DuplicateFilter",
    eventTag = cms.InputTag("eventMaker","evtevent"),
    hypLtP4Tag = cms.InputTag("hypDilepMaker","hypltp4"),
    trksD0Tag = cms.InputTag("trackMaker","trksd0"),
    runTag = cms.InputTag("eventMaker","evtrun"),
    removeEvents = cms.vstring()
)


