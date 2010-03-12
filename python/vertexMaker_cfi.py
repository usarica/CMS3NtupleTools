import FWCore.ParameterSet.Config as cms

vertexMaker = cms.EDFilter("VertexMaker",
	aliasPrefix = cms.untracked.string("vtxs"),
                           primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
)

