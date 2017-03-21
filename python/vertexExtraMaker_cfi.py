import FWCore.ParameterSet.Config as cms

vertexExtraMaker = cms.EDProducer(
  "VertexExtraMaker",
  aliasPrefix           = cms.untracked.string("vtxs"),
  primaryVertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
)

