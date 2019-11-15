import FWCore.ParameterSet.Config as cms

cms3ntuple = cms.EDAnalyzer(
   "CMS3Ntuplizer",

   year = cms.int32(-1), # Must be overriden by main_pset
   treename = cms.untracked.string("Events"),

   isMC = cms.bool(False),

   prefiringWeightsTag = cms.untracked.string(""),

   electronSrc   = cms.InputTag("electronMaker"),
   photonSrc   = cms.InputTag("photonMaker"),
   muonSrc   = cms.InputTag("muonMaker"),
   ak4jetSrc   = cms.InputTag("pfJetMaker"),

   pfmetSrc = cms.InputTag("pfmetMaker"),
   puppimetSrc = cms.InputTag("pfmetpuppiMaker"),

   vtxSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),

   rhoSrc   = cms.InputTag("fixedGridRhoFastjetAll"),
   triggerInfoSrc = cms.InputTag("hltMaker"),
   puInfoSrc = cms.InputTag("slimmedAddPileupInfo"),
   metFilterInfoSrc = cms.InputTag("metFilterMaker"),

   genInfoSrc = cms.InputTag("genMaker"),
   prunedGenParticlesSrc  = cms.InputTag("prunedGenParticles"),
   genJetsSrc  = cms.InputTag("slimmedGenJets"),

   )

