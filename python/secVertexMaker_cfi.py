import FWCore.ParameterSet.Config as cms
import CMS3.NtupleMaker.configProcessName as configProcessName

#from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *

secondaryVertexMaker = cms.EDProducer("SecondaryVertexMaker",

  aliasPrefix             = cms.untracked.string("svs"),
  primaryVertexInputTag   = cms.InputTag("offlineSlimmedPrimaryVertices"),
  inclusiveVertexInputTag = cms.InputTag("slimmedSecondaryVertices","",configProcessName.name),
)

#cms2InclusiveVertexing = cms.Sequence(inclusiveVertexing*secondaryVertexMaker)
