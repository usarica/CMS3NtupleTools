import FWCore.ParameterSet.Config as cms

beamSpotMaker = cms.EDProducer("BeamSpotMaker",
	aliasPrefix = cms.untracked.string("evt_bs"),
    beamSpotInputTag = cms.InputTag("offlineBeamSpot")
)


