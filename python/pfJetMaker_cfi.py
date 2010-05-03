import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDFilter("PFJetMaker",
             pfJetsInputTag     = cms.InputTag("ak5PFJets"),
             pfJetPtCut         = cms.double(5.),
             PFJetCorrectorL2L3 = cms.string("ak5PFL2L3")
)


