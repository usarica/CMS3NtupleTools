import FWCore.ParameterSet.Config as cms

trackToElsAssMaker = cms.EDFilter("TrackToElAssMaker",
    electronsInputTag = cms.InputTag("gsfElectrons"),
    tracksInputTag    = cms.InputTag("generalTracks")
)


