import FWCore.ParameterSet.Config as cms

unpackedTracksAndVertices = cms.EDProducer('PATTrackAndVertexUnpacker',
  slimmedVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  additionalTracks= cms.InputTag("lostTracks"),
  packedCandidates = cms.InputTag("packedPFCandidates"),
  slimmedSecondaryVertices= cms.InputTag("slimmedSecondaryVertices")
)

vertexMaker = cms.EDProducer(
  "VertexMaker",
  aliasPrefix           = cms.untracked.string("vtxs"),
  #primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
  primaryVertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
  #primaryVertexInputTag = cms.InputTag("unpackedTracksAndVertices")
)

vertexMakerWithBS = cms.EDProducer(
  "VertexMaker",
  aliasPrefix           = cms.untracked.string("bsvtxs"),
  primaryVertexInputTag = cms.InputTag("offlinePrimaryVerticesWithBS"),
)
