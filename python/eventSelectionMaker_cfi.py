import FWCore.ParameterSet.Config as cms

eventSelectionMaker = cms.EDProducer("EventSelectionMaker",
	aliasPrefix           = cms.untracked.string("evtsel"),
    primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
    tracksInputTag        = cms.InputTag("generalTracks"),
)

