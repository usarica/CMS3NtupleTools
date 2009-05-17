import FWCore.ParameterSet.Config as cms

jetMaker = cms.EDFilter("JetMaker",
           uncorJetsInputTag = cms.InputTag("sisCone5CaloJets"),
           L2L3corJetsInputTag = cms.InputTag("L2L3CorJet"),
           L2L3L4corJetsInputTag = cms.InputTag("L2L3L4CorJet")
)


