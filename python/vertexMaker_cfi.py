import FWCore.ParameterSet.Config as cms

vertexMaker = cms.EDFilter("VertexMaker",
                           primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
)

