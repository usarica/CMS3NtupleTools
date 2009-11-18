import FWCore.ParameterSet.Config as cms

jetMaker = cms.EDFilter("JetMaker",
           uncorJetsInputTag     = cms.InputTag("prunedUncorrectedCMS2Jets"),
           runningOnReco         = cms.untracked.bool(True),
           L2L3JetCorrectorName = cms.string("L2L3JetCorrectorAK5Calo"),
           AliasPrefix          = cms.string("jets"),
           jetIDIputTag = cms.InputTag("cms2ak5JetID")
)


