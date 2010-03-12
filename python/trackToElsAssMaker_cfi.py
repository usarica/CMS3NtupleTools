import FWCore.ParameterSet.Config as cms

trackToElsAssMaker = cms.EDFilter("TrackToElAssMaker",
	aliasPrefix = cms.untracked.string("trks"),
    electronsInputTag = cms.InputTag("gsfElectrons"),
    tracksInputTag    = cms.InputTag("generalTracks")
)


