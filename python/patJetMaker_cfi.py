import FWCore.ParameterSet.Config as cms

patJetMaker = cms.EDFilter("PATJetMaker",
    # qt jet collection
    patJetsInputTag  = cms.InputTag("selectedPatJets"),
    uncorRecoJetsTag = cms.InputTag("prunedUncorrectedCMS2Jets")
)


