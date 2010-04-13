import FWCore.ParameterSet.Config as cms

scjetMaker = cms.EDFilter("JetMaker",
           uncorJetsInputTag     = cms.InputTag("prunedUncorrectedCMS2scJets"),
           runningOnReco         = cms.untracked.bool(True),
           L2L3JetCorrectorName = cms.string("L2L3JetCorrectorSC5Calo"),
           correctionLevels      = cms.string("L2:L3"),
           correctionTags        = cms.string("Summer09_7TeV_L2Relative_SC5Calo:Summer09_7TeV_L3Absolute_SC5Calo"),
           AliasPrefix          = cms.string("scjets"),
           jetIDIputTag = cms.InputTag("cms2sc5JetID")
)


