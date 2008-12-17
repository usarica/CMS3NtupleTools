import FWCore.ParameterSet.Config as cms

trackToElsAssMaker = cms.EDFilter("TrackToElAssMaker",
    # min DR
    minDR = cms.double(0.1),
    haveHits = cms.bool(True),
    electronsInputTag = cms.InputTag("selectedLayer1Electrons"),
    tracksInputTag    = cms.InputTag("generalTracks")
)


