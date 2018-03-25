import FWCore.ParameterSet.Config as cms

subJetMaker = cms.EDProducer("SubJetMaker",
  # pfJetsInputTag                   = cms.InputTag("prunedUncorrectedCMS2Jets", "pfjet"),
  # pfJetsInputTag                   = cms.InputTag("patJetsAK8PFCHS"),
  # pfJetsInputTag                   = cms.InputTag("ak8PFJetsCHS"),
  pfJetsInputTag                   = cms.InputTag("slimmedJetsAK8"),
  pfCandidatesTag                  = cms.InputTag("packedPFCandidates"),
  pfJetPtCut                       = cms.double(200),
  lessBranches                     = cms.bool(True),

)
