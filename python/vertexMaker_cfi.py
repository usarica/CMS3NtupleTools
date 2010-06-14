import FWCore.ParameterSet.Config as cms

vertexMaker = cms.EDProducer("VertexMaker",
	aliasPrefix = cms.untracked.string("vtxs"),
                           primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
)

