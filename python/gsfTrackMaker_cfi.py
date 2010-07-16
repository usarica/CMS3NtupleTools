import FWCore.ParameterSet.Config as cms

gsfTrackMaker = cms.EDProducer("GSFTrackMaker",
    gsftracksInputTag          = cms.InputTag("electronGsfTracks"),
    beamSpotInputTag        = cms.InputTag("beamSpotMaker", "evtbsp4")
)


