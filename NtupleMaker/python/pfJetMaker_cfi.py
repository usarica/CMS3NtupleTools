import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer(
   "PFJetMaker",

   year = cms.int32(-1), # Must be overriden by main_pset

   aliasprefix = cms.untracked.string("ak4pfjets_chs"),
   jetCollection = cms.untracked.string("AK4PFchs"),

   isMC = cms.bool(False),

   rhoInputTag = cms.InputTag("fixedGridRhoFastjetAll"),
   vtxInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),

   pfJetsInputTag = cms.InputTag("slimmedJets"),
   pfCandidatesInputTag = cms.InputTag("packedPFCandidates"),
   genJetsInputTag = cms.InputTag("slimmedGenJets"),

   METshift_fixEE2017 = cms.bool(False),
   enableManualMETfix = cms.bool(False),
   JEClevels = cms.vstring(),

   )

pfJetPUPPIMaker = cms.EDProducer(
   "PFJetMaker",

   year = cms.int32(-1), # Must be overriden by main_pset

   aliasprefix = cms.untracked.string("ak4pfjets_puppi"),
   jetCollection = cms.untracked.string("AK4PFPuppi"),

   isMC = cms.bool(False),

   rhoInputTag = cms.InputTag("fixedGridRhoFastjetAll"),
   vtxInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),

   pfJetsInputTag = cms.InputTag("slimmedJetsPuppi"),
   pfCandidatesInputTag = cms.InputTag("packedPFCandidates"),
   genJetsInputTag = cms.InputTag("slimmedGenJets"),

   METshift_fixEE2017 = cms.bool(False),
   enableManualMETfix = cms.bool(False),
   JEClevels = cms.vstring(),

   )

subJetMaker = cms.EDProducer(
   "PFJetMaker",

   year = cms.int32(-1), # Must be overriden by main_pset

   aliasprefix = cms.untracked.string("ak8pfjets"),
   jetCollection = cms.untracked.string("AK8PFPuppi"),

   isMC = cms.bool(False),

   rhoInputTag = cms.InputTag("fixedGridRhoFastjetAll"),
   vtxInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),

   pfJetsInputTag = cms.InputTag("slimmedJetsAK8"),
   pfCandidatesInputTag = cms.InputTag("packedPFCandidates"),
   genJetsInputTag = cms.InputTag("slimmedGenJetsAK8"),

   METshift_fixEE2017 = cms.bool(False),
   enableManualMETfix = cms.bool(False),
   JEClevels = cms.vstring(),

   )
