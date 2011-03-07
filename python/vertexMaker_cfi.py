import FWCore.ParameterSet.Config as cms

vertexMaker = cms.EDProducer("VertexMaker",
	aliasPrefix = cms.untracked.string("vtxs"),
                           primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
)

#D.A. vertices will become default at some point.... (4_2_X???)
from CMS2.NtupleMaker.daVtxSequence_cff import *
davertexMaker = cms.EDProducer("VertexMaker",
	aliasPrefix = cms.untracked.string("davtxs"),
                           primaryVertexInputTag = cms.InputTag("offlinePrimaryVerticesDA"),
)

