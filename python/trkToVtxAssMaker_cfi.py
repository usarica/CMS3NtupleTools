import FWCore.ParameterSet.Config as cms

trkToVtxAssMaker = cms.EDFilter("TrkToVtxAssMaker",
	aliasPrefix = cms.untracked.string("trks"),
    vtxInputTag = cms.InputTag("vertexMaker"),
    trksInputTag = cms.InputTag("trackMaker")
)


