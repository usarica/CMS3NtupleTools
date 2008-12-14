import FWCore.ParameterSet.Config as cms

trackMaker = cms.EDFilter("TrackMaker",
    # track collection
    tracksInputTag = cms.InputTag("generalTracks"),
    beamSpotInputTag = cms.InputTag("beamSpotMaker"),
    trkIsolationdRConeMin = cms.double(0.01),
    trkIsolationdRConeMax = cms.double(0.3),
    trkIsolationVtxDiffZMax = cms.double(0.5),
    trkIsolationTkVtxDMax = cms.double(0.1),
    trkIsolationPtMin = cms.double(1.0),
    trkIsolationNHits = cms.int32(7)
)


