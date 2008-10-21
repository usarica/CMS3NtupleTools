import FWCore.ParameterSet.Config as cms

trackMaker = cms.EDFilter("TrackMaker",
    # track collection
    tracksInputTag = cms.InputTag("generalTracks"),
    beamSpotInputTag = cms.InputTag("beamSpotMaker")
)


