import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.DefaultJEC_cff import *

jetMaker = cms.EDProducer("JetMaker",
           uncorJetsInputTag     = cms.InputTag("prunedUncorrectedCMS2Jets", "calojet"),
           AliasPrefix           = cms.string("jets"),
           CaloJetCorrectorL2L3 = cms.string('ak5CaloL2L3')
           runningOnReco         = cms.untracked.bool(True),
           CaloJetCorrectorL1L2L3 = cms.string('ak5CaloL1L2L3'),
           CaloJetCorrectorL1FastL2L3 = cms.string('ak5CaloL1FastL2L3')
)


