import FWCore.ParameterSet.Config as cms

scjetMaker = cms.EDFilter("JetMaker",
           uncorJetsInputTag     = cms.InputTag("prunedUncorrectedCMS2scJets"),
           runningOnReco         = cms.untracked.bool(True),
           L2L3JetCorrectorName = cms.string("L2L3JetCorrectorSC5Calo"),
           AliasPrefix          = cms.string("scjets"),
           jetIDIputTag = cms.InputTag("cms2sc5JetID")
)


