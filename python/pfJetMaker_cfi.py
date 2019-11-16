import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer(
   "PFJetMaker",
   aliasprefix = cms.untracked.string("pfjets"),
   jetCollection = cms.untracked.string("AK4PFchs"),

   isMC = cms.bool(False),

   rhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
   vtxInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),

   pfJetsInputTag = cms.InputTag("slimmedJets"),
   pfCandidatesInputTag = cms.InputTag("packedPFCandidates"),
   genJetsInputTag = cms.InputTag("slimmedGenJets"),

   )

pfJetPUPPIMaker = cms.EDProducer(
   "PFJetMaker",
   aliasprefix = cms.untracked.string("pfjets_puppi"),
   jetCollection = cms.untracked.string("AK4PFPuppi"),

   isMC = cms.bool(False),

   rhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
   vtxInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),

   pfJetsInputTag = cms.InputTag("slimmedJetsPuppi"),
   pfCandidatesInputTag = cms.InputTag("packedPFCandidates"),
   genJetsInputTag = cms.InputTag("slimmedGenJets"),

   )

