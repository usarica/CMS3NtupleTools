import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer("PFJetMaker",
             pfJetsInputTag     = cms.InputTag("prunedUncorrectedCMS2Jets", "pfjet"),
             pfJetPtCut         = cms.double(5.),
             PFJetCorrectorL2L3 = cms.string("ak5PFL2L3")
)


