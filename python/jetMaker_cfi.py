import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.DefaultJEC_cff import *

jetMaker = cms.EDFilter("JetMaker",
           uncorJetsInputTag     = cms.InputTag("prunedUncorrectedCMS2Jets", "calojet"),
           runningOnReco         = cms.untracked.bool(True),
           AliasPrefix          = cms.string("jets"),
           jetIDInputTag = cms.InputTag("cms2ak5JetID"),
           CaloJetCorrectorL2L3 = cms.string('ak5CaloL2L3')
)


