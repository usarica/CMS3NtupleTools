import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer("PFJetMaker",
  aliasPrefix = cms.untracked.string("pfjets"),
  pfJetsInputTag                   = cms.InputTag("slimmedJets"),
  pfJetPtCut                       = cms.double(5.),
  pfCandidatesTag     = cms.InputTag("packedPFCandidates")
)

pfJetPUPPIMaker = cms.EDProducer("PFJetMaker",
  aliasPrefix = cms.untracked.string("pfjets_puppi"),
  pfJetsInputTag                   = cms.InputTag("slimmedJetsPuppi"),
  pfJetPtCut                       = cms.double(20.), #miniAOD doesn't go lower than 20
  pfCandidatesTag     = cms.InputTag("packedPFCandidates")
)
