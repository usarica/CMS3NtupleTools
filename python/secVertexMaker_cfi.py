import FWCore.ParameterSet.Config as cms

from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *

secondaryVertexMaker = cms.EDProducer("SecondaryVertexMaker",
	aliasPrefix = cms.untracked.string("svs"),
        primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
        inclusiveVertexInputTag = cms.InputTag("inclusiveVertices"),
)

cms2InclusiveVertexing = cms.Sequence(inclusiveVertexing*secondaryVertexMaker)
