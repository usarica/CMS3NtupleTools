import FWCore.ParameterSet.Config as cms

trkToVtxAssMaker = cms.EDFilter("TrkToVtxAssMaker",
    vtxInputTag = cms.InputTag("vertexMaker"),
    trksInputTag = cms.InputTag("trackMaker")
)


