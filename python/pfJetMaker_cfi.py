import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer(
   "PFJetMaker",
   aliasprefix = cms.untracked.string("pfjets"),
   jetCollection = cms.untracked.string("AK4PFchs"),

   isMC = cms.bool(False),

   pfJetsInputTag                   = cms.InputTag("slimmedJets"),
   pfJetPtCut                       = cms.double(0.),
   pfCandidatesTag     = cms.InputTag("packedPFCandidates")
   )

pfJetPUPPIMaker = cms.EDProducer(
   "PFJetMaker",
   aliasprefix = cms.untracked.string("pfjets_puppi"),
   jetCollection = cms.untracked.string("AK4PFPuppi"),

   isMC = cms.bool(False),

   pfJetsInputTag                   = cms.InputTag("slimmedJetsPuppi"),
   pfJetPtCut                       = cms.double(0.), #miniAOD doesn't go lower than 20
   pfCandidatesTag     = cms.InputTag("packedPFCandidates")
   )

