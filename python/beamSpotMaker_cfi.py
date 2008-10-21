import FWCore.ParameterSet.Config as cms

beamSpotMaker = cms.EDFilter("BeamSpotMaker",
    beamSpotInputTag = cms.InputTag("offlineBeamSpot")
)


