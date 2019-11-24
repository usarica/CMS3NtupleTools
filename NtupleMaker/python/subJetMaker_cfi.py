import FWCore.ParameterSet.Config as cms

subJetMaker = cms.EDProducer(
   "SubJetMaker",
   # pfJetsInputTag                   = cms.InputTag("prunedUncorrectedCMS2Jets", "pfjet"),
   # pfJetsInputTag                   = cms.InputTag("patJetsAK8PFCHS"),
   # pfJetsInputTag                   = cms.InputTag("ak8PFJetsCHS"),
   jetCollection = cms.untracked.string("AK8PFPuppi"),

   isMC = cms.bool(False),

   pfJetsInputTag = cms.InputTag("slimmedJetsAK8"),
   pfCandidatesTag = cms.InputTag("packedPFCandidates"),
   pfJetPtCut = cms.double(0),
   lessBranches = cms.bool(False),
)
