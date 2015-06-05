import FWCore.ParameterSet.Config as cms

#from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *

secondaryVertexMaker = cms.EDProducer("SecondaryVertexMaker",

  aliasPrefix             = cms.untracked.string("svs"),
  primaryVertexInputTag   = cms.InputTag("offlineSlimmedPrimaryVertices"),
  inclusiveVertexInputTag = cms.InputTag("slimmedSecondaryVertices","","PAT"),
)

#cms2InclusiveVertexing = cms.Sequence(inclusiveVertexing*secondaryVertexMaker)
