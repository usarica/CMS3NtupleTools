import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer("PFJetMaker",
  aliasPrefix = cms.untracked.string("pfjets"),
  pfJetsInputTag                   = cms.InputTag("slimmedJets"),
  pfCandidatesTag                  = cms.InputTag("packedPFCandidates"),
  pfJetPtCut                       = cms.double(5.),
  #PFJetCorrectorL2L3               = cms.string("ak4PFCHSL2L3"),
  #PFJetCorrectorL1Fast             = cms.string("ak4PFCHSL1Fastjet"),
  #PFJetCorrectorL1FastL2L3         = cms.string("ak4PFCHSL1FastL2L3")
  #PFJetCorrectorL1FastL2L3residual = cms.string("ak4PFCHSL1FastL2L3")
)

pfJetPUPPIMaker = cms.EDProducer("PFJetMaker",
  aliasPrefix = cms.untracked.string("pfjets_puppi"),
  pfJetsInputTag                   = cms.InputTag("slimmedJetsPuppi"),
  pfCandidatesTag                  = cms.InputTag("packedPFCandidates"),
  pfJetPtCut                       = cms.double(20.), #miniAOD doesn't go lower than 20
  #PFJetCorrectorL2L3               = cms.string("ak4PFCHSL2L3"),
  #PFJetCorrectorL1Fast             = cms.string("ak4PFCHSL1Fastjet"),
  #PFJetCorrectorL1FastL2L3         = cms.string("ak4PFCHSL1FastL2L3")
  #PFJetCorrectorL1FastL2L3residual = cms.string("ak4PFCHSL1FastL2L3")
)
