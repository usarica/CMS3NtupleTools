import FWCore.ParameterSet.Config as cms

trackMaker = cms.EDProducer("TrackMaker",
	aliasPrefix = cms.untracked.string("trks"),
    tracksInputTag          = cms.InputTag("generalTracks"),
    beamSpotInputTag        = cms.InputTag("beamSpotMaker", "evtbsp4"),
    trkIsolationdRConeMin   = cms.double(0.01),
    trkIsolationdRConeMax   = cms.double(0.3),
    trkIsolationVtxDiffZMax = cms.double(0.5),
    trkIsolationTkVtxDMax   = cms.double(0.1),
    trkIsolationPtMin       = cms.double(1.0),
    trkIsolationNHits       = cms.int32(7)
)


