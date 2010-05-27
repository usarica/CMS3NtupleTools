import FWCore.ParameterSet.Config as cms

patJetMaker = cms.EDFilter("PATJetMaker",
    # qt jet collection
    aliasPrefox = cms.untracked.string("jets_pat"),
    patJetsInputTag  = cms.InputTag("selectedPatJets"),
    uncorRecoJetsTag = cms.InputTag("prunedUncorrectedCMS2Jets", "calojet")
)


