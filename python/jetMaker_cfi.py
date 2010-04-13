import FWCore.ParameterSet.Config as cms

jetMaker = cms.EDFilter("JetMaker",
           uncorJetsInputTag     = cms.InputTag("prunedUncorrectedCMS2Jets"),
           runningOnReco         = cms.untracked.bool(True),
           correctionLevels      = cms.string("L2:L3"),
           correctionTags        = cms.string("Summer09_7TeV_L2Relative_AK5Calo:Summer09_7TeV_L3Absolute_AK5Calo"),
           AliasPrefix          = cms.string("jets"),
           jetIDIputTag = cms.InputTag("cms2ak5JetID")
)


