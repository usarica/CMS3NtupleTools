import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDFilter("PFJetMaker",
             pfJetsInputTag   = cms.InputTag("antikt5PFJets"),
             pfJetPtCut       = cms.double(5.),
             L2L3JetCorrectorName = cms.string("L2L3JetCorrectorAK5PF")
)


