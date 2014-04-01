import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer("PFJetMaker",
  #pfJetsInputTag                   = cms.InputTag("prunedUncorrectedCMS2Jets", "pfjet"),
  pfJetsInputTag                   = cms.InputTag("slimmedJets"),
  #pfCandidatesTag                  = cms.InputTag("particleFlow"),
  pfCandidatesTag                  = cms.InputTag("packedPFCandidates"),
  pfJetPtCut                       = cms.double(5.),
  PFJetCorrectorL2L3               = cms.string("ak5PFCHSL2L3"),
  PFJetCorrectorL1L2L3             = cms.string("ak5PFCHSL1L2L3"),
  PFJetCorrectorL1FastL2L3         = cms.string("ak5PFCHSL1FastL2L3"),
  PFJetCorrectorL1Fast             = cms.string("ak5PFCHSL1Fastjet"),
  PFJetCorrectorL1FastL2L3residual = cms.string("ak5PFCHSL1FastL2L3Residual")
)


