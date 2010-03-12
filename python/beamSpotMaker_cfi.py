import FWCore.ParameterSet.Config as cms

beamSpotMaker = cms.EDFilter("BeamSpotMaker",
	aliasPrefix = cms.untracked.string("evt"),
    beamSpotInputTag = cms.InputTag("offlineBeamSpot")
)


