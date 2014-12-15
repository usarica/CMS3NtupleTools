# import FWCore.ParameterSet.Config as cms

# subJetMaker = cms.EDProducer("SubJetMaker",
#   #pfJetsInputTag                   = cms.InputTag("prunedUncorrectedCMS2Jets", "pfjet"),
#   pfJetsInputTag                   = cms.InputTag("slimmedJetsAK8"),
#   #pfCandidatesTag                  = cms.InputTag("particleFlow"),
#   pfCandidatesTag                  = cms.InputTag("packedPFCandidates"),
#   pfJetPtCut                       = cms.double(5.),
#   PFJetCorrectorL2L3               = cms.string("ak8PFCHSL2L3"),
#   PFJetCorrectorL1L2L3             = cms.string("ak8PFCHSL1L2L3"),
#   PFJetCorrectorL1FastL2L3         = cms.string("ak8PFCHSL1FastL2L3"),
#   PFJetCorrectorL1Fast             = cms.string("ak8PFCHSL1Fastjet"),
#   PFJetCorrectorL1FastL2L3residual = cms.string("ak8PFCHSL1FastL2L3Residual")
# )



import FWCore.ParameterSet.Config as cms

subJetMaker = cms.EDProducer("SubJetMaker",
  #pfJetsInputTag                   = cms.InputTag("prunedUncorrectedCMS2Jets", "pfjet"),
  pfJetsInputTag                   = cms.InputTag("slimmedJetsAK8"),
  # pfJetsInputTag                   = cms.InputTag("patJetsAK8PFCHS"),
  # pfJetsInputTag                   = cms.InputTag("ak8PFJetsCHS"),
  #pfCandidatesTag                  = cms.InputTag("particleFlow"),
  pfCandidatesTag                  = cms.InputTag("packedPFCandidates"),
  pfJetPtCut                       = cms.double(5.),

  # PFJetCorrectorL2L3               = cms.string("ak5PFCHSL2L3"),
  # PFJetCorrectorL1L2L3             = cms.string("ak5PFCHSL1L2L3"),
  # PFJetCorrectorL1FastL2L3         = cms.string("ak5PFCHSL1FastL2L3"),
  # PFJetCorrectorL1Fast             = cms.string("ak5PFCHSL1Fastjet"),
  # PFJetCorrectorL1FastL2L3residual = cms.string("ak5PFCHSL1FastL2L3Residual")
)


