import FWCore.ParameterSet.Config as cms

patJetMaker = cms.EDFilter("PATJetMaker",
    # qt jet collection
    patJetsInputTag  = cms.InputTag("selectedLayer1Jets"),
    uncorRecoJetsTag = cms.InputTag("prunedUncorrectedCMS2Jets")
)


